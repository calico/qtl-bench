#!/usr/bin/env python
"""Shared qtl-bench matching for the single-cell QTL builders.

Every set is built with qtl-bench's own functions (../eqtl_vcfs.py, imported as qb): region
classification, multi-gene collapse, nearest-TSS reassignment, scoring, greedy pick_best, and
exact-allele (SNP) or ilen + mismatch-penalty (indel) buckets. Only the inputs differ from GTEx.

The VCF writers are qtl-bench's, whose INFO descriptions name GTEx; set_info_descriptions()
rewrites those header lines for the single-cell inputs (data lines are untouched).
"""
import os
import sys

import numpy as np
import pandas as pd

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(SRC)
sys.path.insert(0, SRC)
import eqtl_vcfs as qb  # noqa: E402

# the GTEx study's protein-coding GTF (portal.mk), shared so gene sets match across studies
GTF_PC = f"{REPO}/data/gtex11/portal/gencode48_basic_nort_protein.gtf"


def pin_hash_seed():
    """Re-exec under PYTHONHASHSEED=0: set iteration order otherwise varies per run and can flip a
    few near-tied negatives (seen in tenk10k_pbmc/snp/eqtl). Call first in a builder's main()."""
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable] + sys.argv)


def parse_colon_ids(df, indel_t):
    """'chr1:123:A:G' -> chrom/position/ref/alt, qtl-bench's variant filter, MAF."""
    p = df.variant_id.str.split(":", expand=True)
    df = df.assign(chrom=p[0].values, position=p[1].astype(int).values,
                   ref=p[2].values, alt=p[3].values)
    df = qb.filter_variants(df, indel_t)
    df["maf"] = np.minimum(df.af.values, 1.0 - df.af.values)
    return df[df.maf > 0].copy()


def match(pos, neg, bucket_cols, dist_col, allele_weight=1.0):
    """Greedy qtl-bench matching of pos to neg; returns (neg, [(pos row, neg index, score)]).

    dist_col(prow) gives the candidate distance column (None drops the distance term, as for
    eQTL UTR3/CDS).
    """
    neg = neg.reset_index(drop=True)
    groups = {}
    for key, grp in neg.groupby(bucket_cols, sort=False):
        key = key if isinstance(key, tuple) else (key,)
        groups[key] = (grp, qb.pack_events(grp.event.values))

    pos_variants = set(pos.variant_id)
    used, out = set(), []
    for prow in pos.itertuples():
        hit = groups.get(tuple(getattr(prow, c) for c in bucket_cols))
        if hit is None:
            continue
        cand, ev = hit
        dc = dist_col(prow)
        nv, ni, score = qb.pick_best(
            cand, getattr(prow, dc) if dc else np.nan, prow.tpm, prow.maf, dc,
            used, pos_variants, set(), 1.0, 1.0,
            qb.event_mismatch(ev, prow.event), allele_weight)
        if nv is None:
            continue
        used.add(nv)
        out.append((prow, ni, score))
    return neg, out


def build_eqtl_context(pos_raw, cand, gene_tpms, annot, indel_t, out_prefix, label):
    """One cell type, qtl-bench eqtl_vcfs.py main() from 'positives' onward.

    pos_raw  SuSiE rows with pip >= threshold: variant_id, gene, pip, af, afc, pval_nominal
    cand     pool rows already classified and nearest-TSS reassigned (orig_gene kept),
             restricted to this cell type's SuSiE genes, focal positives removed
    """
    utr3_pr, cds_pr, c_utr3_pr, c_cds_pr, gene_tss, _ = annot
    pos = qb.classify_region(pos_raw, utr3_pr, cds_pr, c_utr3_pr, c_cds_pr, gene_tss)
    pos = qb.collapse_multigene(pos).reset_index(drop=True)

    neg = cand.copy()
    neg["tpm"] = neg.gene.map(gene_tpms)
    neg = neg[neg.tpm.notna()]
    neg, pairs = match(pos, neg,
                       ["region", "ref", "alt"] if indel_t == 0 else ["region", "ilen"],
                       lambda p: "tss_dist" if p.region == "TSS" else None)

    rows = []
    for p, ni, score in pairs:
        n = neg.loc[ni]
        rows.append(dict(
            pos_variant=p.variant_id, neg_variant=n.variant_id, pos_gene=p.gene,
            neg_gene=n.gene, region=p.region, pos_tss_dist=p.tss_dist,
            neg_tss_dist=n.tss_dist, pos_tpm=p.tpm, neg_tpm=n.tpm, pos_maf=p.maf,
            neg_maf=n.maf, pos_pip=p.pip,
            mismatch=int(qb.event_mismatch(qb.pack_events([n.event]), p.event)[0]),
            score=score))
    m = pd.DataFrame(rows)
    m.to_csv(f"{out_prefix}_matches.tsv", sep="\t", index=False)
    qb.write_vcf(f"{out_prefix}_pos.vcf", pos, label)
    nw = neg.loc[[ni for _, ni, _ in pairs]].assign(pip=0.0)
    qb.write_vcf(f"{out_prefix}_neg.vcf", nw, label)
    return pos, m


def set_info_descriptions(vcf_file, desc):
    """Rewrite the Description of the ##INFO lines whose ID is in desc; data lines untouched."""
    with open(vcf_file) as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if not line.startswith("#"):
            break
        if line.startswith("##INFO=<ID="):
            iid = line[len("##INFO=<ID="):].split(",", 1)[0]
            if iid in desc:
                lines[i] = line[:line.index('Description="')] + f'Description="{desc[iid]}">\n'
    with open(vcf_file, "w") as fh:
        fh.writelines(lines)


def merge_vcfs(out_file, files):
    """qtl-bench make_vcfs.merge_variants, sorted in pandas (bedtools is not on PATH here)."""
    files = sorted(files)
    with open(files[0]) as fh:
        header = [line for line in fh if line.startswith("#")]
    rows, seen = [], set()
    for f in files:
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                vid = line.split("\t", 3)[2]
                if vid not in seen:
                    seen.add(vid)
                    c = line.split("\t", 2)
                    rows.append((c[0], int(c[1]), line))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_file, "w") as out:
        out.writelines(header)
        out.writelines(r[2] for r in rows)
    return len(rows)


def merge_matches(out_file, files, suffix="_matches.tsv"):
    dfs = []
    for f in sorted(files):
        d = pd.read_csv(f, sep="\t")
        d.insert(0, "tissue", os.path.basename(f)[:-len(suffix)])
        dfs.append(d)
    if dfs:
        pd.concat(dfs, ignore_index=True).to_csv(out_file, sep="\t", index=False)
