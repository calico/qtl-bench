#!/usr/bin/env python
import argparse
import os
import time

import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm

"""
paqtl_vcfs.py

Generate positive and negative apaQTL sets from GTEx v11 SuSiE fine-mapping.
Positives: SuSiE_summary with pip >= pos_pip, protein-coding genes.
  PAS distance = distance to nearest transcript 3' end of the associated gene
  (gene-specific, from gencode GTF transcript records).
Negatives: signif_pairs with high-PIP (variant, gene) pairs removed, matched
  to positives on PAS distance, expression TPM, MAF, and alleles (exact for
  SNPs; same indel type and length, with a per-base mismatch penalty).
  PAS distance = genome-wide nearest transcript 3' end; associated gene is
  reassigned to the gene owning that nearest site (used for TPM lookup).
"""


################################################################################
# main
################################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Generate positive/negative apaQTL sets for one GTEx tissue."
    )
    parser.add_argument("susie_parquet", help="GTEx apaQTLs.SuSiE_summary parquet file")
    parser.add_argument(
        "--signif", default=None,
        help="signif_pairs parquet (default: inferred from susie_parquet)")
    parser.add_argument(
        "-g", dest="genes_gtf",
        default="data/gtex11/portal/gencode48_basic_nort_protein.gtf",
        help="Protein-coding GTF for PAS (transcript 3' end) extraction")
    parser.add_argument(
        "--gtex_expr",
        default="data/gtex11/portal/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct",
        help="GTEx median TPM GCT file")
    parser.add_argument("-o", dest="out_prefix", default="gtex_paqtl")
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
    parser.add_argument("--signif_pool", default=None,
                        help="Cross-tissue pooled signif parquet; replaces --signif as the negative candidate source")
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
    susie_df["gene"] = susie_df.gene_id.str.split(".").str[0]
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

    pc_genes = set(susie_df.gene.unique())

    ############################################################################
    # build PAS annotation (transcript 3' ends) from GTF
    t0 = time.time()
    print(f"[{tissue}] Building PAS annotation from GTF...", flush=True)
    gene_pas_sites, all_pas_pr = build_pas_sites(args.genes_gtf)
    n_total = sum(len(v) for v in gene_pas_sites.values())
    print(f"  {len(gene_pas_sites):,} genes, {n_total:,} unique sites in {time.time()-t0:.0f}s")

    ############################################################################
    # load signif_pairs (candidate negative pool)
    t0 = time.time()
    pool_path = args.signif_pool or args.signif
    print(f"[{tissue}] Loading signif_pairs from {pool_path}...", flush=True)
    pool_cols = ["phenotype_id", "group_id", "variant_id", "af"]
    if not using_pool:
        pool_cols += ["pval_nominal", "slope"]
    signif_df = pd.read_parquet(pool_path, columns=pool_cols)
    if using_pool:
        signif_df["pval_nominal"] = np.nan
        signif_df["slope"] = np.nan
    signif_df["gene"] = signif_df.group_id.str.split(".").str[0]
    signif_df = signif_df[signif_df.gene.isin(pc_genes)].copy()
    signif_df = parse_variant_id(signif_df)
    signif_df = filter_variants(signif_df, args.indel_t)
    signif_df["maf"] = np.minimum(signif_df.af.values, 1.0 - signif_df.af.values)
    signif_df["tpm"] = signif_df.gene.map(gene_tpms)
    signif_df = signif_df[signif_df.tpm.notna()].copy()
    print(f"  {len(signif_df):,} rows after biotype+tpm filter in {time.time()-t0:.0f}s")

    # compute genome-wide PAS distances for unique variants in signif_df
    t0 = time.time()
    print(f"[{tissue}] Computing PAS distances for signif_pairs...", flush=True)
    uniq_vdf = signif_df.drop_duplicates("variant_id")[["variant_id", "chrom", "position"]].copy()
    pas_dists, pas_genes = compute_pas_dist_negatives(uniq_vdf, all_pas_pr)
    pas_dist_map = dict(zip(uniq_vdf.variant_id, pas_dists))
    pas_gene_map = dict(zip(uniq_vdf.variant_id, pas_genes))
    signif_df["pas_dist"] = signif_df.variant_id.map(pas_dist_map).astype(float)
    print(f"  PAS distances computed in {time.time()-t0:.0f}s")

    ############################################################################
    # positives: pip >= pos_pip, gene-specific PAS distances
    pos_raw = susie_df[susie_df.pip >= args.pos_pip].copy()
    print(f"[{tissue}] {len(pos_raw)} positive (variant, phenotype) rows, "
          f"{pos_raw.variant_id.nunique()} unique variants")

    pos_raw["pas_dist"] = compute_pas_dist_positives(pos_raw, gene_pas_sites)
    n_no_site = int(pos_raw.pas_dist.isna().sum())
    if n_no_site:
        print(f"  dropped {n_no_site} positive rows with no PAS for gene")
        pos_raw = pos_raw[pos_raw.pas_dist.notna()].copy()

    # one row per variant: keep maximum-pip row
    pos_df = (
        pos_raw.sort_values("pip", ascending=False, kind="stable")
        .drop_duplicates("variant_id")
        .copy()
    )
    print(f"  {len(pos_df)} positive variants after collapsing")

    # pval and slope lookup for positives — always from per-tissue signif_pairs
    pos_vars_list = list(set(pos_df.variant_id))
    pval_lkp = pd.read_parquet(
        args.signif,
        columns=["phenotype_id", "variant_id", "pval_nominal", "slope"],
        filters=[("variant_id", "in", pos_vars_list)],
    )
    pval_map = dict(zip(zip(pval_lkp.phenotype_id, pval_lkp.variant_id), pval_lkp.pval_nominal))
    slope_map = dict(zip(zip(pval_lkp.phenotype_id, pval_lkp.variant_id), pval_lkp.slope))
    pos_df["pval_nominal"] = [pval_map.get((pid, vid))
                              for pid, vid in zip(pos_df.phenotype_id, pos_df.variant_id)]
    pos_df["slope"] = [slope_map.get((pid, vid))
                       for pid, vid in zip(pos_df.phenotype_id, pos_df.variant_id)]
    del pval_lkp

    pos_variants = set(pos_df.variant_id)

    exclude_variants = set()
    if args.exclude is not None:
        with open(args.exclude) as f:
            exclude_variants = set(f.read().split())
        print(f"  loaded {len(exclude_variants)} excluded variants")

    ############################################################################
    # filter candidate negatives
    # cross-tissue blocklist subsumes focal-tissue positives; filter explicitly
    # for safety when --exclude is not provided (standalone use).
    neg_df = signif_df[~signif_df.variant_id.isin(pos_variants)].copy()
    if exclude_variants:
        neg_df = neg_df[~neg_df.variant_id.isin(exclude_variants)].copy()

    # reassign gene to nearest PAS gene and update TPM for negatives
    neg_df = neg_df.copy()
    neg_df["gene"] = neg_df.variant_id.map(pas_gene_map)
    neg_df["tpm"] = neg_df.gene.map(gene_tpms)
    neg_df = neg_df[neg_df.gene.notna() & neg_df.tpm.notna()].copy()

    # deduplicate to one row per variant: keep minimum pas_dist
    neg_df = (
        neg_df.sort_values("pas_dist")
        .drop_duplicates("variant_id")
        .reset_index(drop=True)
    )
    print(f"  {len(neg_df):,} candidate negatives after filtering")

    ############################################################################
    # matching loop
    # SNPs bucket on exact alleles (only 12 combinations, so exactness is free).
    # Indel alleles are far too diverse for that, so they bucket on signed length --
    # fixing the event type and size, which is what displaces the model's output
    # bins -- and pay an allele_weight penalty per mismatched base in the score.
    bucket_cols = ["ref", "alt"] if args.indel_t == 0 else ["ilen"]
    neg_groups = {
        key: (grp, pack_events(grp.event.values))
        for key, grp in neg_df.groupby(bucket_cols, sort=False)
    }

    neg_variants = set()
    matches = []
    for prow in tqdm(list(pos_df.itertuples()), desc=f"[{tissue}] Matching"):
        key = (prow.ref, prow.alt) if args.indel_t == 0 else (prow.ilen,)
        if key not in neg_groups:
            continue
        cand, cand_events = neg_groups[key]
        nv, ni, score = pick_best(
            cand, prow.pas_dist, prow.tpm, prow.maf,
            "pas_dist", neg_variants, pos_variants, exclude_variants,
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
                "pos_pas_dist": prow.pas_dist,
                "neg_pas_dist": neg_row.pas_dist,
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

    write_vcf(f"{args.out_prefix}_pos.vcf", pos_df, pos_variants, "pas_dist", "PD", tissue)
    neg_write_df = neg_df[neg_df.variant_id.isin(neg_variants)].copy()
    neg_write_df["pip"] = 0.0
    write_vcf(f"{args.out_prefix}_neg.vcf", neg_write_df, neg_variants, "pas_dist", "PD", tissue)


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
# helpers: PAS annotation
################################################################################
def build_pas_sites(genes_gtf):
    """Extract poly-A site positions (transcript 3' ends) from a gencode GTF.

    For + strand transcripts: PAS = transcript end (GTF 1-based inclusive).
    For - strand transcripts: PAS = transcript start (GTF 1-based inclusive).

    Returns:
        gene_pas_sites: dict gene_id (versionless) -> np.int64 array of 1-based positions
        all_pas_pr: pr.PyRanges with gene_id column for genome-wide k_nearest
    """
    gtf = pd.read_csv(
        genes_gtf, sep="\t", comment="#", header=None,
        names=["seqname", "source", "feature", "start", "end",
               "score", "strand", "frame", "attributes"],
        dtype={"seqname": str},
    )
    tx = gtf[gtf.feature == "transcript"].copy()
    tx["gene_id"] = (
        tx.attributes.str.extract(r'gene_id "([^"]+)"')[0].str.split(".").str[0]
    )
    tx["tx_start"] = tx.start.astype(int)
    tx["tx_end"] = tx.end.astype(int)
    # 3' end: GTF is 1-based inclusive; end for + strand, start for - strand
    tx["pas_pos"] = np.where(tx.strand == "+", tx.tx_end, tx.tx_start)

    site_df = (
        tx[["seqname", "pas_pos", "gene_id"]]
        .rename(columns={"seqname": "chrom", "pas_pos": "pos"})
        .drop_duplicates(["chrom", "pos", "gene_id"])
        .copy()
    )
    site_df["pos"] = site_df.pos.astype(int)

    gene_pas_sites = {
        gid: grp.pos.values.astype(np.int64)
        for gid, grp in site_df.groupby("gene_id")
    }

    all_pas_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": site_df.chrom,
        "Start": site_df.pos.astype(int) - 1,  # 0-based half-open
        "End": site_df.pos.astype(int),
        "gene_id": site_df.gene_id,
    }))

    return gene_pas_sites, all_pas_pr


def compute_pas_dist_positives(pos_df, gene_pas_sites):
    """Distance to nearest transcript 3' end of the variant's associated gene.

    Returns a list of floats (np.nan when gene has no entries in gene_pas_sites).
    """
    distances = []
    for row in pos_df.itertuples():
        sites = gene_pas_sites.get(row.gene)
        if sites is None or len(sites) == 0:
            distances.append(np.nan)
        else:
            distances.append(float(np.abs(sites - row.position).min()))
    return distances


def compute_pas_dist_negatives(variants_df, all_pas_pr):
    """Genome-wide nearest PAS (transcript 3' end) distance for a set of unique variants.

    Uses pyranges k_nearest. Assigns gene_id from the nearest PAS.

    Returns (distances, gene_ids) aligned to variants_df row order.
    """
    v_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": variants_df.chrom.values,
        "Start": variants_df.position.values - 1,
        "End": variants_df.position.values,
        "Name": variants_df.variant_id.values,
    }))
    result_df = v_pr.k_nearest(all_pas_pr, ties="first", apply_strand_suffix=False).as_df()
    result_df["abs_dist"] = result_df.Distance.abs()
    best = (
        result_df.sort_values("abs_dist")
        .drop_duplicates("Name")
        .set_index("Name")
    )
    distances = []
    gene_ids = []
    for v in variants_df.variant_id:
        if v in best.index:
            distances.append(float(best.loc[v, "abs_dist"]))
            gene_ids.append(best.loc[v, "gene_id"].split(".")[0])
        else:
            distances.append(np.nan)
            gene_ids.append(None)
    return distances, gene_ids


################################################################################
# helpers: scoring
################################################################################
def compute_scores(candidate_df, target_dist, target_tpm, target_maf,
                   dist_col, dist_weight, maf_weight,
                   mismatch=None, allele_weight=1.0):
    """Score = dist_weight * log1p(dist) L1 + log1p(TPM) L1 + maf_weight * log(MAF) L1.

    Lower is better. MAF term omitted when target_maf is NaN.
    """
    log_dist_diff = np.abs(
        np.log1p(candidate_df[dist_col].values) - np.log1p(target_dist))
    log_tpm_diff = np.abs(
        np.log1p(candidate_df["tpm"].values) - np.log1p(target_tpm))
    scores = dist_weight * log_dist_diff + log_tpm_diff
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
def write_vcf(vcf_file, df, variants_set, dist_col, dist_tag, tissue):
    """Write a VCF with apaQTL-specific INFO fields."""
    header = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=GENE,Number=1,Type=String,Description="apaQTL gene (ENSG, versionless)">',
        '##INFO=<ID=TISSUE,Number=1,Type=String,Description="GTEx tissue">',
        '##INFO=<ID=PIP,Number=1,Type=Float,Description="SuSiE posterior inclusion probability (0 for negatives)">',
        f'##INFO=<ID={dist_tag},Number=1,Type=Integer,Description="PAS distance (nearest 3-prime end of associated gene for positives; genome-wide nearest for negatives)">',
        '##INFO=<ID=TPM,Number=1,Type=Float,Description="GTEx median TPM of gene in tissue">',
        '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">',
        '##INFO=<ID=NLP,Number=1,Type=Float,Description="Negative log10 nominal p-value from tensorQTL">',
        '##INFO=<ID=SLOPE,Number=1,Type=Float,Description="Per-allele change in poly-A site usage fraction from tensorQTL">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]),
    ]
    df = df[df.variant_id.isin(variants_set)].copy()
    if "pip" in df.columns:
        df = df.sort_values("pip", ascending=False, kind="stable")

    with open(vcf_file, "w") as out:
        out.write("\n".join(header) + "\n")
        seen = set()
        for v in df.itertuples():
            if v.variant_id in seen:
                continue
            seen.add(v.variant_id)
            pip = float(getattr(v, "pip", 0.0))
            dist_val = getattr(v, dist_col)
            dist_str = "." if (dist_val is None or (isinstance(dist_val, float) and np.isnan(dist_val))) else str(int(dist_val))
            pval = getattr(v, "pval_nominal", None)
            nlp = "." if (pval is None or (isinstance(pval, float) and np.isnan(pval))) else f"{min(-np.log10(max(pval, 1e-300)), 300):.4f}"
            slope_val = getattr(v, "slope", None)
            slope_str = "." if (slope_val is None or (isinstance(slope_val, float) and np.isnan(slope_val))) else f"{slope_val:.4f}"
            info = (
                f"GENE={v.gene};TISSUE={tissue};PIP={pip:.4f};"
                f"{dist_tag}={dist_str};TPM={v.tpm:.4f};MAF={v.maf:.4f};"
                f"NLP={nlp};SLOPE={slope_str}"
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
