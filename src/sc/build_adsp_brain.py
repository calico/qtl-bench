#!/usr/bin/env python
"""adsp_brain: matched brain snuc-eQTL sets from ADSP FunGen (NG00184), qtl-bench recipe.

Supersedes build_brain_xqtl_sets.py, which reimplemented the matching (TSS-only regions, no
allele buckets, credible-set-member pool) because it assumed FunGen had no signif_pairs. It
does: bh_significant is the per-(variant, gene) BH-significant table. Inputs now mirror GTEx:

  positives   single_context_finemapping_all (scfmAll) PIP >= 0.9, protein-coding target
              (gene in the qtl-bench GTF), AFC = conditional_effect (as before), p = qtl_p_val
  pool        bh_significant pairs of all six cell types in the cohort (qtl-bench
              --signif_pool), af = median non_ref_af, restricted to the focal cell type's
              scfmAll genes (qtl-bench pc_genes)
  blocklist   any variant with PIP >= 0.01 for a protein-coding gene in any cell type of any
              of the three cohorts (qtl-bench high_pip_variants.txt; cohorts share donors)
  expression  --expr siletti (release): HCA Human Brain Cell Atlas v1.0 (Siletti et al. 2023)
              pseudobulk CPM over the atlas clusters of each cell type (--siletti_cpm, genes x
              cell types); --expr mssm: microglia-atlas (MSSM PD multiome) pseudobulk CPM, both
              neuron types on its one Neurons column (--expr_dir, one {column}.tsv per class)
  annotation  GENCODE 48 basic protein-coding GTF, as qtl-bench
  matching    qb_sets: qtl-bench functions, plain GTEx matching
  cohort      one set per ROSMAP cohort (--cohort): CUIMC1+2+MIT -> eqtl/, MIT -> eqtl_rosmap_mit/,
              CUIMC1 -> eqtl_rosmap_cuimc1/; all three feed the blocklist (they share donors)

  python src/sc/build_adsp_brain.py --indel_t {0,4} --cohort C --siletti_cpm CPM   # -> data/adsp_brain/{snp,indel}/<set>/
"""
import argparse
import glob
import json
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

import qb_sets as qs
from qb_sets import qb

PORTAL = f"{qs.REPO}/data/adsp_brain/portal"  # fetch_adsp_brain.sh puts raw/NG00184/ here
COHORTS = ["ROSMAP_MIT", "ROSMAP_CUIMC1", "ROSMAP_CUIMC1_2_MIT"]
SET_NAME = {"ROSMAP_CUIMC1_2_MIT": "eqtl", "ROSMAP_MIT": "eqtl_rosmap_mit", "ROSMAP_CUIMC1": "eqtl_rosmap_cuimc1"}
# FunGen context code -> (cell type name used downstream, microglia atlas column)
CT = {
    "ast": ("astrocyte", "Astrocytes"),
    "exc": ("excitatory_neuron", "Neurons"),
    "inh": ("inhibitory_neuron", "Neurons"),
    "mic": ("microglia", "Myeloid"),
    "oli": ("oligodendrocyte", "Oligodendrocytes"),
    "opc": ("OPC", "OPCs"),
}
POS_PIP, EXCLUDE_PIP = 0.9, 0.01
INFO_DESC = {
    "TISSUE": "ADSP FunGen xQTL snuc-eQTL cell type (ROSMAP)",
    "AFC": "FunGen SuSiE conditional_effect (alt vs ref); positives only",
    "NLP": "Negative log10 nominal p-value (FunGen qtl_p_val)",
}
TPM_DESC = {
    "mssm": "MSSM microglia-atlas pseudobulk CPM of gene; neuron types share the Neurons column",
    "siletti": "HCA brain atlas v1.0 pseudobulk CPM of gene over the cell type's clusters",
}
FM_COLS = ["variant_id", "PIP", "target_ensembl_id", "conditional_effect", "non_ref_af",
           "target_info"]
BH_COLS = ["variant_id", "target_ensembl_id", "non_ref_af"]


def files(raw, cohort, kind, ct):
    sub = {"fm": "single_context_finemapping_all", "bh": "bh_significant"}[kind]
    return sorted(glob.glob(f"{raw}/{cohort}/snuc-eQTL/{sub}/bed3plus17_qtl/hg38/"
                            f"*.snuc-eQTL.{ct}.*.bed.gz"))


def read_fm(path):
    d = pd.read_csv(path, sep="\t", usecols=FM_COLS)
    d["pval_nominal"] = pd.to_numeric(
        d.target_info.str.extract(r"qtl_p_val=([^;]+)")[0], errors="coerce")
    d["non_ref_af"] = pd.to_numeric(d.non_ref_af, errors="coerce")   # "NR" -> NaN, dropped
    return d.drop(columns="target_info")


def read_bh(path):
    d = pd.read_csv(path, sep="\t", usecols=BH_COLS)
    d["non_ref_af"] = pd.to_numeric(d.non_ref_af, errors="coerce")
    return d


def load(raw, kind, cohort, ct, workers):
    with Pool(workers) as p:
        parts = p.map(read_fm if kind == "fm" else read_bh, files(raw, cohort, kind, ct))
    d = pd.concat(parts, ignore_index=True)
    return d.rename(columns={"target_ensembl_id": "gene", "non_ref_af": "af",
                             "PIP": "pip", "conditional_effect": "afc"})


def main():
    qs.pin_hash_seed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--indel_t", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cohort", default="ROSMAP_CUIMC1_2_MIT", choices=COHORTS)
    ap.add_argument("--expr", choices=["mssm", "siletti"], default="siletti")
    ap.add_argument("--expr_dir", help="--expr mssm: directory of {atlas column}.tsv (gene_id, expression)")
    ap.add_argument("--siletti_cpm", help="--expr siletti: pseudobulk CPM table, genes x cell types")
    ap.add_argument("--portal", default=PORTAL, help="directory holding raw/NG00184/")
    ap.add_argument("--out", help="output set directory (default data/adsp_brain/{snp,indel}/eqtl)")
    ap.add_argument("--gtf", default=qs.GTF_PC, help="qtl-bench protein-coding GTF")
    args = ap.parse_args()
    np.random.seed(0)
    vclass = "snp" if args.indel_t == 0 else "indel"
    raw_dir = f"{args.portal}/raw/NG00184/ADSP_FunGen_xQTL/v1"
    out = args.out or f"{qs.REPO}/data/adsp_brain/{vclass}/{SET_NAME[args.cohort]}"

    print("Region annotation...", flush=True)
    annot = qb.build_region_annotation(args.gtf)
    utr3_pr, cds_pr, c_utr3_pr, c_cds_pr, gene_tss, all_tss_pr = annot
    known = set(gene_tss) | set(utr3_pr.df.gene_id) | set(cds_pr.df.gene_id)
    if args.expr == "mssm":  # keyed by benchmark cell type name
        cols = {col: dict(zip(d.gene_id.str.split(".").str[0], d.expression))
                for col in {c for _, c in CT.values()}
                for d in [pd.read_csv(f"{args.expr_dir}/{col}.tsv", sep="\t")]}
        expr = {name: cols[col] for name, col in CT.values()}
    else:
        cpm = pd.read_csv(args.siletti_cpm, sep="\t", index_col=0)
        expr = {name: cpm[name].to_dict() for name, _ in CT.values()}
    expr_genes = set().union(*[set(v) for v in expr.values()])

    # fine-mapping for every (cohort, cell type): positives and the cross-cohort blocklist
    fm = {}
    for cohort in COHORTS:
        for ct in CT:
            d = load(raw_dir, "fm", cohort, ct, args.workers)
            d["gene"] = d.gene.str.split(".").str[0]
            fm[cohort, ct] = d
            print(f"  scfmAll {cohort} {ct}: {len(d):,} rows", flush=True)
    block = set().union(*[set(d.loc[(d.pip >= EXCLUDE_PIP) & d.gene.isin(known), "variant_id"])
                          for d in fm.values()])
    print(f"  blocklist: {len(block):,} variants with PIP >= {EXCLUDE_PIP}", flush=True)

    summary = []
    for cohort in [args.cohort]:
        os.makedirs(out, exist_ok=True)

        # pooled bh_significant for the cohort; classification and nearest-TSS reassignment
        # depend only on (variant, gene), so they run once
        pool = pd.concat([load(raw_dir, "bh", cohort, ct, args.workers) for ct in CT])
        n_raw = len(pool)
        pool["gene"] = pool.gene.str.split(".").str[0]
        pool = pool[pool.gene.isin(known) & ~pool.variant_id.isin(block)]
        pool = pool.groupby(["variant_id", "gene"], as_index=False, sort=False).af.median()
        pool = qs.parse_colon_ids(pool, args.indel_t)
        print(f"[{cohort}] pool {n_raw:,} bh rows -> {len(pool):,} (variant, gene) candidates",
              flush=True)
        pool = qb.classify_region(pool, utr3_pr, cds_pr, c_utr3_pr, c_cds_pr, gene_tss)
        pool["orig_gene"] = pool.gene
        pool["tpm"] = 0.0
        # placeholder TPM map: drop TSS negatives whose nearest gene has no atlas row
        pool = qb.reassign_tss_negatives(pool, all_tss_pr, {g: 0.0 for g in expr_genes})
        pool = pool.reset_index(drop=True)
        print(f"  regions {pool.region.value_counts().to_dict()}", flush=True)
        info = dict(INFO_DESC, TPM=TPM_DESC[args.expr])

        for ct, (name, col) in CT.items():
            gene_tpms = expr[name]
            sus = fm[cohort, ct]
            sus = sus[sus.gene.isin(known)]
            sus = sus.assign(tpm=sus.gene.map(gene_tpms))
            sus = sus[sus.tpm.notna()]
            pc_genes = set(sus.gene)
            raw = fm[cohort, ct]
            pos_raw = sus[sus.pip >= POS_PIP]
            pos_raw = qs.parse_colon_ids(pos_raw, args.indel_t)
            cand = pool[pool.orig_gene.isin(pc_genes) & ~pool.variant_id.isin(set(pos_raw.variant_id))]
            pos, m = qs.build_eqtl_context(pos_raw, cand, gene_tpms, annot, args.indel_t,
                                           f"{out}/{name}", name)
            for side in ("pos", "neg"):
                qs.set_info_descriptions(f"{out}/{name}_{side}.vcf", info)
            row = dict(
                cohort=cohort, cell_type=name,
                pip90_rows=int((raw.pip >= POS_PIP).sum()),
                pip90_variants=raw.loc[raw.pip >= POS_PIP, "variant_id"].nunique(),
                protein_coding_expressed=sus.loc[sus.pip >= POS_PIP, "variant_id"].nunique(),
                **{vclass: len(pos)}, matched=len(m),
                candidates=len(cand), median_score=float(m.score.median()) if len(m) else np.nan)
            summary.append(row)
            print(f"  {name:<18} {row}", flush=True)

        qs.merge_vcfs(f"{out}/merge.vcf",
                      glob.glob(f"{out}/*_pos.vcf") + glob.glob(f"{out}/*_neg.vcf"))
        qs.merge_matches(f"{out}/matches.tsv", glob.glob(f"{out}/*_matches.tsv"))
        with open(f"{out}/build.json", "w") as fh:
            json.dump({"pip": POS_PIP, "exclude_pip": EXCLUDE_PIP, "indel_t": args.indel_t,
                       "positives": "scfmAll PIP >= pip, protein-coding, atlas-expressed",
                       "pool": "bh_significant, pooled over 6 cell types, focal SuSiE genes",
                       "blocklist": "PIP >= exclude_pip, any cohort / cell type",
                       "expression": {"mssm": "microglia atlas (MSSM) pseudobulk CPM",
                                      "siletti": "HCA brain atlas v1.0 pseudobulk CPM, track clusters"}[args.expr],
                       "gtf": os.path.basename(args.gtf),
                       "per_cell_type": [r for r in summary if r["cohort"] == cohort]},
                      fh, indent=2)


if __name__ == "__main__":
    main()
