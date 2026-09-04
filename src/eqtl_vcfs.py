#!/usr/bin/env python
import argparse
import os
import time

import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm

"""
eqtl_vcfs.py

Generate positive and negative eQTL sets from GTEx v11 SuSiE fine-mapping.
Positives come from SuSiE_summary (pip >= pos_pip, protein-coding). Negatives
come from signif_pairs, with any (variant, gene) pair whose SuSiE PIP meets
exclude_pip removed, matched to positives within a region bucket (UTR3, CDS,
or TSS-distance) on expression TPM, MAF, distance to TSS, and alleles
(exact for SNPs; same indel type and length, with a per-base mismatch penalty).
"""

REGION_RANK = {"UTR3": 0, "CDS": 1, "TSS": 2}


################################################################################
# main
################################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Generate positive/negative eQTL sets for one GTEx tissue."
    )
    parser.add_argument("susie_parquet", help="GTEx SuSiE_summary parquet file")
    parser.add_argument(
        "--signif", default=None,
        help="signif_pairs parquet (default: inferred from susie_parquet)")
    parser.add_argument(
        "--signif_pool", default=None,
        help="Pre-pooled signif_pairs parquet (variant_id, phenotype_id, af) across all "
             "tissues; when set, replaces --signif as the negative candidate source.")
    parser.add_argument(
        "-g", dest="genes_gtf",
        default="data/gtex11/portal/gencode48_basic_nort_protein.gtf",
        help="Protein-coding GTF for UTR/CDS/TSS intervals")
    parser.add_argument(
        "--gtex_expr",
        default="data/gtex11/portal/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct",
        help="GTEx median TPM GCT file")
    parser.add_argument("-o", dest="out_prefix", default="gtex_eqtl")
    parser.add_argument("--pos_pip", default=0.9, type=float,
                        help="PIP threshold for positives")
    parser.add_argument("--exclude_pip", default=0.01, type=float,
                        help="PIP above which a (variant, gene) pair is removed from the negative pool")
    parser.add_argument("--indel_t", default=0, type=int,
                        help="0 = SNPs only; >0 = indels only, up to this length")
    parser.add_argument("--dist_weight", default=1.0, type=float)
    parser.add_argument("--maf_weight", default=1.0, type=float)
    parser.add_argument("--allele_weight", default=1.0, type=float,
                        help="Score penalty per mismatched base between the "
                             "positive and negative indel alleles (indels only)")
    parser.add_argument("--exclude", default=None,
                        help="File of variant IDs to never use as negatives")
    args = parser.parse_args()

    np.random.seed(0)

    tissue = os.path.basename(args.susie_parquet).split(".v11.")[0]
    if args.signif is None:
        args.signif = args.susie_parquet.replace("SuSiE_summary", "signif_pairs")
    using_pool = args.signif_pool is not None

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    ############################################################################
    # load SuSiE summary, filter to protein-coding, parse variant_id
    t0 = time.time()
    print(f"[{tissue}] Loading SuSiE summary...", flush=True)
    susie_df = pd.read_parquet(args.susie_parquet)
    susie_df["gene"] = susie_df.phenotype_id.str.split(".").str[0]
    susie_df = susie_df[susie_df.biotype == "protein_coding"].copy()
    susie_df = parse_variant_id(susie_df)
    susie_df = filter_variants(susie_df, args.indel_t)
    susie_df["maf"] = np.minimum(susie_df.af.values, 1.0 - susie_df.af.values)
    print(f"  {len(susie_df):,} protein-coding SuSiE rows in {time.time()-t0:.0f}s")

    ############################################################################
    # GTEx TPM
    gene_tpms = load_tpms(args.gtex_expr, tissue)
    susie_df["tpm"] = susie_df.gene.map(gene_tpms)
    susie_df = susie_df[susie_df.tpm.notna()].copy()

    ############################################################################
    # region annotation (shared between positives and candidates)
    t0 = time.time()
    print(f"[{tissue}] Building region annotation...", flush=True)
    utr3_pr, cds_pr, canonical_utr3_pr, canonical_cds_pr, gene_tss, all_tss_pr = build_region_annotation(args.genes_gtf)
    print(f"  DONE in {time.time()-t0:.0f}s")

    # drop genes absent from annotation (e.g. protein_coding_LoF-only genes
    # filtered out of gencode48_basic_nort_protein.gtf); their TSS distance would be
    # NaN and poison the matching scores.
    known_genes = set(gene_tss.keys())
    known_genes |= set(utr3_pr.df.gene_id.unique()) if len(utr3_pr) else set()
    known_genes |= set(cds_pr.df.gene_id.unique()) if len(cds_pr) else set()
    before_rows = len(susie_df)
    before_genes = susie_df.gene.nunique()
    susie_df = susie_df[susie_df.gene.isin(known_genes)].copy()
    dropped_rows = before_rows - len(susie_df)
    dropped_genes = before_genes - susie_df.gene.nunique()
    if dropped_rows:
        print(f"  dropped {dropped_rows:,} rows ({dropped_genes} genes) absent from annotation")

    pc_genes = set(susie_df.gene.unique())

    ############################################################################
    # positives
    pos_raw = susie_df[susie_df.pip >= args.pos_pip].copy()
    # drop rows where GTEx SuSiE produced no AFC estimate (rare; concentrated
    # in segmental-duplication / paralog-cluster genes where SuSiE signals
    # collapse across multiple credible sets and AFC is not emitted).
    afc_na = pos_raw.afc.isna()
    if afc_na.any():
        print(f"  dropping {afc_na.sum()} positive rows with no AFC")
        pos_raw = pos_raw[~afc_na].copy()
    print(f"[{tissue}] {len(pos_raw)} positive (variant, gene) rows, "
          f"{pos_raw.variant_id.nunique()} unique variants")
    pos_raw = classify_region(pos_raw, utr3_pr, cds_pr, canonical_utr3_pr, canonical_cds_pr, gene_tss)
    pos_df = collapse_multigene(pos_raw)
    print(f"  after multi-gene collapse: {len(pos_df)} positive variants")
    print(f"  region counts: {pos_df.region.value_counts().to_dict()}")

    # pval lookup for positives — always from per-tissue signif_pairs
    pos_vars_list = list(set(pos_df.variant_id))
    pval_lkp = pd.read_parquet(
        args.signif,
        columns=["phenotype_id", "variant_id", "pval_nominal"],
        filters=[("variant_id", "in", pos_vars_list)],
    )
    pval_map = dict(zip(zip(pval_lkp.phenotype_id, pval_lkp.variant_id), pval_lkp.pval_nominal))
    pos_df["pval_nominal"] = [pval_map.get((pid, vid))
                              for pid, vid in zip(pos_df.phenotype_id, pos_df.variant_id)]
    del pval_lkp

    pos_variants = set(pos_df.variant_id)

    exclude_variants = set()
    if args.exclude is not None:
        with open(args.exclude) as f:
            exclude_variants = set(f.read().split())
        print(f"  loaded {len(exclude_variants)} excluded variants")

    ############################################################################
    # candidate negative pool from signif_pairs (per-tissue) or signif_pool (cross-tissue)
    t0 = time.time()
    pool_path = args.signif_pool or args.signif
    pool_label = "signif_pool" if args.signif_pool else "signif_pairs"
    print(f"[{tissue}] Loading {pool_label}...", flush=True)
    pool_cols = ["phenotype_id", "variant_id", "af"]
    if not using_pool:
        pool_cols.append("pval_nominal")
    signif_df = pd.read_parquet(pool_path, columns=pool_cols)
    if using_pool:
        signif_df["pval_nominal"] = np.nan
    signif_df["gene"] = signif_df.phenotype_id.str.split(".").str[0]
    # gene must be expressed in this tissue and have a SuSiE row (defines pc_genes)
    signif_df = signif_df[signif_df.gene.isin(pc_genes)].copy()
    signif_df = parse_variant_id(signif_df)
    signif_df = filter_variants(signif_df, args.indel_t)
    signif_df["maf"] = np.minimum(signif_df.af.values, 1.0 - signif_df.af.values)
    signif_df["tpm"] = signif_df.gene.map(gene_tpms)
    signif_df = signif_df[signif_df.tpm.notna()].copy()
    print(f"  {len(signif_df):,} rows after biotype+tpm filter")

    # cross-tissue blocklist subsumes focal-tissue positives (PIP ≥ 0.9 ≥ exclude_pip),
    # but filter explicitly for safety when --exclude is not provided (standalone use).
    signif_df = signif_df[~signif_df.variant_id.isin(pos_variants)].copy()
    if exclude_variants:
        signif_df = signif_df[~signif_df.variant_id.isin(exclude_variants)].copy()
    print(f"  {len(signif_df):,} candidates after filtering, "
          f"loaded in {time.time()-t0:.0f}s")

    ############################################################################
    # classify candidate regions
    t0 = time.time()
    print(f"[{tissue}] Classifying candidate regions...", flush=True)
    neg_df = classify_region(signif_df, utr3_pr, cds_pr, canonical_utr3_pr, canonical_cds_pr, gene_tss)
    print(f"  DONE in {time.time()-t0:.0f}s")
    print(f"  region counts: {neg_df.region.value_counts().to_dict()}")

    # For TSS-region negatives: replace gene/tss_dist with genome-wide nearest TSS.
    # The signif_pairs-associated gene may have a distant TSS while the variant is
    # actually much closer to a different gene's TSS.
    t0 = time.time()
    print(f"[{tissue}] Reassigning TSS distances for TSS-region negatives...", flush=True)
    neg_df = reassign_tss_negatives(neg_df, all_tss_pr, gene_tpms)
    print(f"  DONE in {time.time()-t0:.0f}s, "
          f"{(neg_df.region == 'TSS').sum()} TSS-region candidates remaining")

    ############################################################################
    # matching loop
    neg_df = neg_df.reset_index(drop=True)
    neg_variants = set()
    matches = []

    # Pre-group candidates so each positive's lookup is O(|group|). SNPs bucket on
    # exact alleles (only 12 combinations, so exactness is free). Indel alleles are
    # far too diverse for that, so they bucket on signed length -- fixing the event
    # type and size, which is what displaces the model's output bins -- and pay an
    # allele_weight penalty per mismatched base in the score instead.
    bucket_cols = ["region", "ref", "alt"] if args.indel_t == 0 else ["region", "ilen"]
    neg_groups = {
        key: (grp, pack_events(grp.event.values))
        for key, grp in neg_df.groupby(bucket_cols, sort=False)
    }
    if neg_groups:
        sizes = np.array([len(g) for g, _ in neg_groups.values()])
        print(f"  {len(neg_groups)} {tuple(bucket_cols)} buckets; "
              f"size median={int(np.median(sizes))} p10={int(np.quantile(sizes, 0.1))} "
              f"max={int(sizes.max())}")

    for prow in tqdm(list(pos_df.itertuples()), desc=f"[{tissue}] Matching"):
        key = ((prow.region, prow.ref, prow.alt) if args.indel_t == 0
               else (prow.region, prow.ilen))
        if key not in neg_groups:
            continue
        cand, cand_events = neg_groups[key]
        dist_col = "tss_dist" if prow.region == "TSS" else None
        nv, ni, score = pick_best(
            cand, prow.tss_dist, prow.tpm, prow.maf,
            dist_col, neg_variants, pos_variants, exclude_variants,
            args.dist_weight, args.maf_weight,
            event_mismatch(cand_events, prow.event), args.allele_weight,
        )
        if nv is not None:
            neg_variants.add(nv)
            neg_row = neg_df.loc[ni]
            matches.append({
                "pos_variant": prow.variant_id,
                "neg_variant": nv,
                "pos_gene": prow.gene,
                "neg_gene": neg_row.gene,
                "region": prow.region,
                "pos_tss_dist": prow.tss_dist,
                "neg_tss_dist": neg_row.tss_dist,
                "pos_tpm": prow.tpm,
                "neg_tpm": neg_row.tpm,
                "pos_maf": prow.maf,
                "neg_maf": neg_row.maf,
                "pos_pip": prow.pip,
                "mismatch": int(event_mismatch(
                    pack_events([neg_row.event]), prow.event)[0]),
                "score": score,
            })

    print(f"[{tissue}] {len(pos_variants)} positives, {len(neg_variants)} matched negatives")
    if len(neg_variants) < len(pos_variants):
        print(f"  WARNING: {len(pos_variants) - len(neg_variants)} unmatched positives")

    matches_df = pd.DataFrame(matches)
    matches_df.to_csv(f"{args.out_prefix}_matches.tsv", sep="\t", index=False)

    write_vcf(f"{args.out_prefix}_pos.vcf", pos_df, tissue)
    neg_write_df = neg_df[neg_df.variant_id.isin(neg_variants)].copy()
    neg_write_df["pip"] = 0.0
    # keep the (variant, gene) row that was actually selected
    selected = matches_df[["neg_variant", "neg_gene"]].rename(
        columns={"neg_variant": "variant_id", "neg_gene": "gene"})
    neg_write_df = neg_write_df.merge(selected, on=["variant_id", "gene"], how="inner")
    neg_write_df = neg_write_df.drop_duplicates("variant_id")
    write_vcf(f"{args.out_prefix}_neg.vcf", neg_write_df, tissue)


################################################################################
# helpers: input parsing
################################################################################
def parse_variant_id(df):
    """Parse GTEx variant_id 'chr1_64764_C_T_b38' into chrom/position/ref/alt."""
    parts = df.variant_id.str.split("_", expand=True)
    df = df.copy()
    df["chrom"] = parts[0].values
    df["position"] = parts[1].astype(int).values
    df["ref"] = parts[2].values
    df["alt"] = parts[3].values
    return df


def filter_variants(df, indel_t):
    """indel_t == 0: SNPs only. indel_t > 0: indels only, length <= indel_t.

    Adds two derived columns used by matching: `ilen`, the signed length change
    (>0 insertion, <0 deletion, 0 SNP), and `event`, the longer of (ref, alt) --
    the inserted or deleted bases plus their anchor base.
    """
    lr = df.ref.str.len()
    la = df.alt.str.len()
    if indel_t == 0:
        df = df[(lr == 1) & (la == 1)].copy()
    else:
        df = df[(lr - la).abs().between(1, indel_t)].copy()
    df["ilen"] = df.alt.str.len() - df.ref.str.len()
    df["event"] = np.where(df.ilen <= 0, df.ref, df.alt)
    return df


def load_tpms(gtex_gct, tissue_col):
    """Load GTEx median TPM for a single tissue column into a gene_id -> TPM dict."""
    gct_df = pd.read_csv(gtex_gct, sep="\t", skiprows=2, usecols=["Name", tissue_col])
    gct_df["gene_id"] = gct_df["Name"].str.split(".").str[0]
    return dict(zip(gct_df.gene_id, gct_df[tissue_col].values))


################################################################################
# helpers: region annotation
################################################################################
def build_region_annotation(genes_gtf):
    """Parse GTF for CDS, 3' UTR, and per-transcript TSS positions.

    Returns:
        utr3_pr:          pyranges of 3' UTR intervals (all basic transcripts) with 'gene_id'
        cds_pr:           pyranges of CDS intervals (all basic transcripts) with 'gene_id'
        canonical_utr3_pr: pyranges of 3' UTR intervals for canonical transcript per gene
        canonical_cds_pr:  pyranges of CDS intervals for canonical transcript per gene
        gene_tss: dict gene_id -> np.ndarray of TSS positions (1-based)
    """
    gtf = pd.read_csv(
        genes_gtf,
        sep="\t", comment="#", header=None,
        names=["seqname", "source", "feature", "start", "end",
               "score", "strand", "frame", "attributes"],
        dtype={"seqname": str},
    )
    gtf = gtf[gtf.feature.isin(["CDS", "UTR", "transcript"])].copy()
    gtf["gene_id"] = gtf.attributes.str.extract(r'gene_id "([^"]+)"')[0].str.split(".").str[0]
    gtf["transcript_id"] = gtf.attributes.str.extract(r'transcript_id "([^"]+)"')[0]

    # Canonical transcript per gene for UTR3/CDS disambiguation:
    # MANE_Select preferred, GENCODE_Primary fallback for genes without it.
    gtf["is_mane"] = gtf.attributes.str.contains('tag "MANE_Select"', regex=False)
    gtf["is_primary"] = gtf.attributes.str.contains('tag "GENCODE_Primary"', regex=False)
    tx = gtf[gtf.feature == "transcript"]
    mane_tx = tx[tx.is_mane][["gene_id", "transcript_id"]]
    primary_tx = tx[tx.is_primary][["gene_id", "transcript_id"]]
    covered = set(mane_tx.gene_id)
    canonical_ids = set(mane_tx.transcript_id) | set(
        primary_tx[~primary_tx.gene_id.isin(covered)].transcript_id)

    cds_rows = gtf[gtf.feature == "CDS"]
    cds_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": cds_rows.seqname.values,
        "Start": cds_rows.start.values.astype(int) - 1,
        "End": cds_rows.end.values.astype(int),
        "gene_id": cds_rows.gene_id.values,
    }))

    # per-transcript CDS span for 3'UTR classification
    cds_span = cds_rows.groupby("transcript_id").agg(
        cds_min=("start", "min"), cds_max=("end", "max"))

    utr_rows = gtf[gtf.feature == "UTR"].merge(
        cds_span, left_on="transcript_id", right_index=True, how="left")
    utr_rows = utr_rows[utr_rows.cds_min.notna()]
    is_utr3_plus = (utr_rows.strand == "+") & (utr_rows.start > utr_rows.cds_max)
    is_utr3_minus = (utr_rows.strand == "-") & (utr_rows.end < utr_rows.cds_min)
    utr3 = utr_rows[is_utr3_plus | is_utr3_minus]
    utr3_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": utr3.seqname.values,
        "Start": utr3.start.values.astype(int) - 1,
        "End": utr3.end.values.astype(int),
        "gene_id": utr3.gene_id.values,
    }))

    # Canonical-transcript-only PyRanges for UTR3/CDS disambiguation
    canonical_cds_rows = cds_rows[cds_rows.transcript_id.isin(canonical_ids)]
    canonical_cds_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": canonical_cds_rows.seqname.values,
        "Start": canonical_cds_rows.start.values.astype(int) - 1,
        "End": canonical_cds_rows.end.values.astype(int),
        "gene_id": canonical_cds_rows.gene_id.values,
    }))

    canonical_cds_span = canonical_cds_rows.groupby("transcript_id").agg(
        cds_min=("start", "min"), cds_max=("end", "max"))
    canonical_utr_rows = gtf[
        (gtf.feature == "UTR") & gtf.transcript_id.isin(canonical_ids)
    ].merge(canonical_cds_span, left_on="transcript_id", right_index=True, how="left")
    canonical_utr_rows = canonical_utr_rows[canonical_utr_rows.cds_min.notna()]
    is_cu3_plus  = (canonical_utr_rows.strand == "+") & (canonical_utr_rows.start > canonical_utr_rows.cds_max)
    is_cu3_minus = (canonical_utr_rows.strand == "-") & (canonical_utr_rows.end   < canonical_utr_rows.cds_min)
    canonical_utr3 = canonical_utr_rows[is_cu3_plus | is_cu3_minus]
    canonical_utr3_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": canonical_utr3.seqname.values,
        "Start": canonical_utr3.start.values.astype(int) - 1,
        "End": canonical_utr3.end.values.astype(int),
        "gene_id": canonical_utr3.gene_id.values,
    }))

    # per-transcript TSS: 1-based position of the first transcribed base, matching
    # the convention of gencode_tss.py -d 1 -u 1 (TSS = transcript.start - 1 for
    # + strand, transcript.end - 1 for - strand).
    tx_rows = gtf[gtf.feature == "transcript"]
    tss_pos = np.where(
        tx_rows.strand.values == "+",
        tx_rows.start.values.astype(np.int64) - 1,
        tx_rows.end.values.astype(np.int64) - 1,
    )
    tx_tss = pd.DataFrame({"gene_id": tx_rows.gene_id.values, "tss": tss_pos})
    gene_tss = {
        gid: grp.tss.values.astype(np.int64)
        for gid, grp in tx_tss.groupby("gene_id")
    }

    # flat PyRanges of all TSS sites for genome-wide k_nearest on negatives
    all_tss_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": tx_rows.seqname.values,
        "Start": tss_pos,
        "End": tss_pos + 1,
        "gene_id": tx_rows.gene_id.values,
    }))

    return utr3_pr, cds_pr, canonical_utr3_pr, canonical_cds_pr, gene_tss, all_tss_pr


def reassign_tss_negatives(neg_df, all_tss_pr, gene_tpms):
    """For TSS-region negatives, replace gene/tss_dist with genome-wide nearest TSS.

    A negative's signif_pairs gene may have a distant TSS even when the variant is
    near a different gene's TSS. Using the genome-wide nearest TSS mirrors how
    sQTL/paQTL negatives are handled and avoids inflated distance estimates.
    Gene is reassigned to the owner of the nearest TSS for TPM lookup.
    """
    tss_mask = neg_df.region.values == "TSS"
    if not tss_mask.any():
        return neg_df

    neg_df = neg_df.copy()
    uniq_vdf = (
        neg_df[tss_mask]
        .drop_duplicates("variant_id")[["variant_id", "chrom", "position"]]
        .copy()
    )

    v_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": uniq_vdf.chrom.values,
        "Start": uniq_vdf.position.values - 1,
        "End": uniq_vdf.position.values,
        "Name": uniq_vdf.variant_id.values,
    }))
    result_df = v_pr.k_nearest(all_tss_pr, ties="first", apply_strand_suffix=False).as_df()
    result_df["abs_dist"] = result_df.Distance.abs()
    best = (
        result_df.sort_values("abs_dist")
        .drop_duplicates("Name")
        .set_index("Name")
    )

    dist_map = {
        v: float(best.loc[v, "abs_dist"]) if v in best.index else np.nan
        for v in uniq_vdf.variant_id
    }
    gene_map = {
        v: best.loc[v, "gene_id"].split(".")[0] if v in best.index else None
        for v in uniq_vdf.variant_id
    }

    tss_idx = neg_df.index[tss_mask]
    neg_df.loc[tss_idx, "tss_dist"] = neg_df.loc[tss_idx, "variant_id"].map(dist_map).values
    neg_df.loc[tss_idx, "gene"] = neg_df.loc[tss_idx, "variant_id"].map(gene_map).values
    neg_df.loc[tss_idx, "tpm"] = neg_df.loc[tss_idx, "gene"].map(gene_tpms).values

    drop_mask = tss_mask & neg_df.tpm.isna()
    if drop_mask.any():
        neg_df = neg_df[~drop_mask].copy()

    return neg_df


def classify_region(df, utr3_pr, cds_pr, canonical_utr3_pr, canonical_cds_pr, gene_tss):
    """Annotate df with 'region' and 'tss_dist' columns.

    Region is relative to the gene in the row.  The union of all basic-annotation
    transcripts determines whether a variant is in a functional category (UTR3 or
    CDS) vs TSS.  When a variant overlaps both UTR3 and CDS across different
    transcripts, the canonical transcript (MANE_Select / GENCODE_Primary) breaks
    the tie; if the canonical transcript also covers neither, UTR3 takes precedence.
    Rows reset to sequential index.
    """
    df = df.reset_index(drop=True).copy()
    n = len(df)

    # span the full ref allele so multi-base deletions are placed correctly;
    # identical to a single base for SNPs
    v_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": df.chrom.values,
        "Start": df.position.values - 1,
        "End": df.position.values - 1 + df.ref.str.len().values,
        "_idx": np.arange(n),
        "gene_id": df.gene.values,
    }))

    any_utr3 = _gene_matched_overlap_idx(v_pr, utr3_pr)
    any_cds  = _gene_matched_overlap_idx(v_pr, cds_pr)

    ambiguous    = any_utr3 & any_cds
    utr3_unambig = any_utr3 - ambiguous
    cds_unambig  = any_cds  - ambiguous

    # Resolve ambiguous variants using the canonical transcript
    canon_utr3    = _gene_matched_overlap_idx(v_pr, canonical_utr3_pr) & ambiguous
    canon_cds     = (_gene_matched_overlap_idx(v_pr, canonical_cds_pr) & ambiguous) - canon_utr3
    canon_neither = ambiguous - canon_utr3 - canon_cds  # fall back to UTR3 precedence

    utr3_set = utr3_unambig | canon_utr3 | canon_neither
    cds_set  = cds_unambig  | canon_cds

    region = np.full(n, "TSS", dtype=object)
    if utr3_set:
        region[np.fromiter(utr3_set, dtype=np.int64)] = "UTR3"
    if cds_set:
        region[np.fromiter(cds_set, dtype=np.int64)] = "CDS"
    df["region"] = region

    # Vectorized TSS distance for TSS rows, grouped by gene
    tss_dist = np.full(n, np.nan)
    tss_mask = df.region.values == "TSS"
    if tss_mask.any():
        tss_df = df[tss_mask]
        positions = tss_df.position.values
        idx = tss_df.index.values
        for gene_id, grp_idx in tss_df.groupby("gene").indices.items():
            tsss = gene_tss.get(gene_id)
            if tsss is None or len(tsss) == 0:
                continue
            sub_pos = positions[grp_idx]
            dists = np.abs(sub_pos[:, None] - tsss[None, :]).min(axis=1)
            tss_dist[idx[grp_idx]] = dists
    df["tss_dist"] = tss_dist
    return df


def _gene_matched_overlap_idx(v_pr, region_pr):
    """Return the set of _idx values whose variant overlaps a region of the same gene_id."""
    if len(region_pr) == 0:
        return set()
    joined = v_pr.join(region_pr, suffix="_b").as_df()
    if len(joined) == 0:
        return set()
    matched = joined[joined.gene_id == joined.gene_id_b]
    return set(matched._idx.values.tolist())


def collapse_multigene(pos_raw):
    """Collapse multi-gene positive rows to one row per variant.

    Precedence: UTR3 > CDS > TSS (ties broken by smallest tss_dist).
    """
    pos_raw = pos_raw.copy()
    pos_raw["_rank"] = pos_raw.region.map(REGION_RANK)
    pos_raw["_tss_sort"] = pos_raw.tss_dist.fillna(1e12)
    pos_raw = pos_raw.sort_values(["variant_id", "_rank", "_tss_sort"])
    return pos_raw.drop_duplicates("variant_id", keep="first").drop(
        columns=["_rank", "_tss_sort"])


################################################################################
# helpers: scoring
################################################################################
def compute_scores(candidate_df, target_dist, target_tpm, target_maf,
                   dist_col, dist_weight, maf_weight,
                   mismatch=None, allele_weight=1.0):
    """Score = log1p(TPM) L1 + [dist_weight * log1p(dist) L1] + [maf_weight * log(MAF) L1].

    Lower is better. Distance term omitted when dist_col is None (UTR3/CDS).
    MAF term omitted when target_maf is NaN.
    """
    log_tpm_diff = np.abs(
        np.log1p(candidate_df["tpm"].values) - np.log1p(target_tpm))
    scores = log_tpm_diff
    if dist_col is not None:
        log_dist_diff = np.abs(
            np.log1p(candidate_df[dist_col].values) - np.log1p(target_dist))
        scores = scores + dist_weight * log_dist_diff
    if not np.isnan(target_maf) and maf_weight > 0 and "maf" in candidate_df.columns:
        log_maf_diff = np.abs(
            np.log(candidate_df["maf"].values) - np.log(target_maf))
        scores = scores + maf_weight * log_maf_diff
    if mismatch is not None:
        scores = scores + allele_weight * mismatch
    return scores


def pack_events(events):
    """Pack equal-length event alleles into an (n, L) byte matrix."""
    return np.frombuffer("".join(events).encode(), dtype="S1").reshape(len(events), -1)


def event_mismatch(packed, target_event):
    """Hamming distance from target_event for each packed candidate event allele.

    Candidates share a bucket's ilen and hence its event length, so this is well
    defined. It is zero whenever the alleles are identical, which for SNP buckets
    (keyed on exact ref/alt) is always.
    """
    return (packed != np.frombuffer(target_event.encode(), dtype="S1")).sum(axis=1)


def pick_best(candidate_df, target_dist, target_tpm, target_maf,
              dist_col, neg_variants, pos_variants, exclude_variants,
              dist_weight, maf_weight, mismatch=None, allele_weight=1.0):
    """Return (variant_id, df_index, score) for the best available candidate."""
    if len(candidate_df) == 0:
        return None, None, None
    scores = compute_scores(
        candidate_df, target_dist, target_tpm, target_maf,
        dist_col, dist_weight, maf_weight, mismatch, allele_weight)
    variant_ids = candidate_df.variant_id.values
    indices = candidate_df.index.values
    for ni in np.argsort(scores, kind="stable"):
        vid = variant_ids[ni]
        if (vid not in neg_variants
                and vid not in pos_variants
                and vid not in exclude_variants):
            return vid, indices[ni], scores[ni]
    return None, None, None


################################################################################
# VCF output
################################################################################
def write_vcf(vcf_file, df, tissue):
    """Write a VCF with eQTL-specific INFO fields."""
    header = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=GENE,Number=1,Type=String,Description="eQTL gene (ENSG, versionless)">',
        '##INFO=<ID=TISSUE,Number=1,Type=String,Description="GTEx tissue">',
        '##INFO=<ID=PIP,Number=1,Type=Float,Description="SuSiE posterior inclusion probability (0 for negatives)">',
        '##INFO=<ID=REGION,Number=1,Type=String,Description="Variant region: UTR3, CDS, or TSS">',
        '##INFO=<ID=TSSD,Number=1,Type=Integer,Description="Distance to nearest TSS of eQTL gene">',
        '##INFO=<ID=TPM,Number=1,Type=Float,Description="GTEx median TPM of gene in tissue">',
        '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">',
        '##INFO=<ID=AFC,Number=1,Type=Float,Description="GTEx SuSiE allelic fold change (log2, alt vs ref); positives only">',
        '##INFO=<ID=NLP,Number=1,Type=Float,Description="Negative log10 nominal p-value from tensorQTL">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]),
    ]
    with open(vcf_file, "w") as out:
        out.write("\n".join(header) + "\n")
        if "pip" in df.columns:
            df = df.sort_values("pip", ascending=False, kind="stable")
        seen = set()
        for v in df.itertuples():
            if v.variant_id in seen:
                continue
            seen.add(v.variant_id)
            pip = float(getattr(v, "pip", 0.0))
            tssd = "." if (v.tss_dist is None or np.isnan(v.tss_dist)) else str(int(v.tss_dist))
            afc_val = getattr(v, "afc", None)
            afc = "." if (afc_val is None or (isinstance(afc_val, float) and np.isnan(afc_val))) else f"{afc_val:.4f}"
            pval = getattr(v, "pval_nominal", None)
            nlp = "." if (pval is None or (isinstance(pval, float) and np.isnan(pval))) else f"{min(-np.log10(max(pval, 1e-300)), 300):.4f}"
            info = (
                f"GENE={v.gene};TISSUE={tissue};PIP={pip:.4f};"
                f"REGION={v.region};TSSD={tssd};"
                f"TPM={v.tpm:.4f};MAF={v.maf:.4f};AFC={afc};NLP={nlp}"
            )
            out.write("\t".join([
                v.chrom, str(v.position), v.variant_id,
                v.ref, v.alt, ".", ".", info,
            ]) + "\n")


################################################################################
# __main__
################################################################################
if __name__ == "__main__":
    main()
