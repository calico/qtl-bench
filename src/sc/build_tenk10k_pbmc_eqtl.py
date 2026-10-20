#!/usr/bin/env python
"""tenk10k_pbmc eqtl: matched TenK10K PBMC sc-eQTL sets, qtl-bench recipe (plain GTEx matching).

SNP (--indel_t 0) and indel (--indel_t 4) sets; negatives are drawn through qb_sets.match.

  positives   TenK10K SuSiE PIP >= 0.9 (finemapped_eqtls_pbmc_celltypes_TenK10K_sc.tsv),
              protein-coding target in the qtl-bench GTF, expressed in the --expr atlas; AF, p and BETA from
              nominal/{ct}_pos_beta.tsv (BETA written to AFC)
  pool        per-gene BH q < 0.05 pairs (nominal/{ct}_bh_pool.tsv.gz), pooled across every built
              cell type (qtl-bench --signif_pool), restricted to focal SuSiE genes
  cell types  --expr onek1k: all 28 TenK10K cell types (each has a OneK1K Azimuth l2 profile);
              --expr hpa: the 7 HPA profiles
  blocklist   any variant with PIP >= 0.01 in any of the 28 cell types
  expression  --expr onek1k (release): OneK1K pseudobulk CPM on the Azimuth l2 labels TenK10K
              uses (--onek1k_cpm, genes x cell types); --expr hpa: HPA rna_immune_cell nTPM (--hpa)
  annotation  GENCODE 48 basic protein-coding GTF
  inputs      --portal holds finemapped_eqtls_pbmc_celltypes_TenK10K_sc.tsv
              (build_pbmc_finemapped_tsv.py) and nominal/ (extract_pbmc_eqtl_nominal.py)

  python src/sc/build_tenk10k_pbmc_eqtl.py --indel_t 0 --onek1k_cpm CPM   # -> data/tenk10k_pbmc/snp/eqtl/
  python src/sc/build_tenk10k_pbmc_eqtl.py --indel_t 4 --onek1k_cpm CPM   # -> data/tenk10k_pbmc/indel/eqtl/
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

import qb_sets as qs
from qb_sets import qb

PORTAL = f"{qs.REPO}/data/tenk10k_pbmc/portal/eqtl"
HPA_CT = {
    "CD4_Naive": "naive CD4 T-cell",
    "NK": "NK-cell",
    "CD4_TCM": "memory CD4 T-cell",
    "CD8_TEM": "memory CD8 T-cell",
    "B_naive": "naive B-cell",
    "CD14_Mono": "classical monocyte",
    "CD16_Mono": "non-classical monocyte",
}
HPA_CT_INV = {v: k for k, v in HPA_CT.items()}
POS_PIP, EXCLUDE_PIP = 0.9, 0.01
INFO_DESC = {
    "TISSUE": "TenK10K cell type (Azimuth l2)",
    "AFC": "TenK10K SAIGE-QTL BETA (alt allele); positives only",
    "NLP": "Negative log10 SAIGE-QTL nominal p-value",
}
TPM_DESC = {
    "hpa": "HPA rna_immune_cell nTPM of gene in the matching immune cell",
    "onek1k": "OneK1K pseudobulk CPM of gene (Azimuth l2)",
}


def main():
    qs.pin_hash_seed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--indel_t", type=int, default=0)
    ap.add_argument("--expr", choices=["hpa", "onek1k"], default="onek1k")
    ap.add_argument("--hpa", help="--expr hpa: HPA rna_immune_cell.tsv")
    ap.add_argument("--onek1k_cpm", help="--expr onek1k: pseudobulk CPM table, genes x cell types")
    ap.add_argument("--portal", default=PORTAL, help="TenK10K eQTL inputs (fine-mapping tsv, nominal/)")
    ap.add_argument("--out", help="output set directory (default data/tenk10k_pbmc/{snp,indel}/eqtl)")
    ap.add_argument("--gtf", default=qs.GTF_PC, help="qtl-bench protein-coding GTF")
    args = ap.parse_args()
    np.random.seed(0)
    vclass = "snp" if args.indel_t == 0 else "indel"
    E = args.portal
    out = args.out or f"{qs.REPO}/data/tenk10k_pbmc/{vclass}/eqtl"
    os.makedirs(out, exist_ok=True)
    info = dict(INFO_DESC, TPM=TPM_DESC[args.expr])

    if args.expr == "hpa":
        hpa = pd.read_csv(args.hpa, sep="\t")
        tpms = {HPA_CT_INV[c]: dict(zip(g.Gene, g.nTPM)) for c, g in hpa.groupby("Immune cell") if c in HPA_CT_INV}
        cts = list(HPA_CT)  # the cell types HPA profiles
        expr_genes = set(hpa.Gene)
    else:
        cpm = pd.read_csv(args.onek1k_cpm, sep="\t", index_col=0)
        col = {c.replace(" ", "_"): c for c in cpm.columns}  # "CD4 Naive" -> CD4_Naive; NK_CD56bright as is
        fm_cts = pd.read_csv(f"{E}/finemapped_eqtls_pbmc_celltypes_TenK10K_sc.tsv", sep="\t", usecols=["cellType"]).cellType
        cts = sorted(c for c in fm_cts.unique() if c in col)  # every TenK10K cell type OneK1K profiles
        tpms = {ct: cpm[col[ct]].to_dict() for ct in cts}
        expr_genes = set(cpm.index)

    print("Region annotation...", flush=True)
    annot = qb.build_region_annotation(args.gtf)
    utr3_pr, cds_pr, c_utr3_pr, c_cds_pr, gene_tss, all_tss_pr = annot
    known = set(gene_tss) | set(utr3_pr.df.gene_id) | set(cds_pr.df.gene_id)

    fm = pd.read_csv(f"{E}/finemapped_eqtls_pbmc_celltypes_TenK10K_sc.tsv", sep="\t",
                     usecols=["variantId", "targetEnsemblID", "cellType", "pip"])
    block = set(fm.loc[fm.pip >= EXCLUDE_PIP, "variantId"])
    print(f"  blocklist: {len(block):,} variants with PIP >= {EXCLUDE_PIP}", flush=True)

    parts = []
    for ct in cts:
        d = pd.read_csv(f"{E}/nominal/{ct}_bh_pool.tsv.gz", sep="\t",
                        usecols=["variant", "gene", "AF_Allele2"])
        parts.append(d[d.gene.isin(known)])
    pool = pd.concat(parts).rename(columns={"variant": "variant_id", "AF_Allele2": "af"})
    n_raw = len(pool)
    pool = pool[~pool.variant_id.isin(block)]
    pool = pool.groupby(["variant_id", "gene"], as_index=False, sort=False).af.median()
    pool = qs.parse_colon_ids(pool, args.indel_t)
    print(f"  pool {n_raw:,} rows -> {len(pool):,} (variant, gene) candidates", flush=True)
    pool = qb.classify_region(pool, utr3_pr, cds_pr, c_utr3_pr, c_cds_pr, gene_tss)
    pool["orig_gene"] = pool.gene
    pool["tpm"] = 0.0
    pool = qb.reassign_tss_negatives(pool, all_tss_pr, {g: 0.0 for g in expr_genes})
    pool = pool.reset_index(drop=True)
    print(f"  regions {pool.region.value_counts().to_dict()}", flush=True)

    summary = []
    for ct in cts:
        gene_tpms = tpms[ct]
        raw = fm[fm.cellType == ct]
        sus = raw[raw.targetEnsemblID.isin(known)].rename(
            columns={"variantId": "variant_id", "targetEnsemblID": "gene"})
        sus = sus.assign(tpm=sus.gene.map(gene_tpms))
        sus = sus[sus.tpm.notna()]
        pc_genes = set(sus.gene)

        nom = pd.read_csv(f"{E}/nominal/{ct}_pos_beta.tsv", sep="\t").rename(
            columns={"variantId": "variant_id", "targetEnsemblID": "gene", "AF_Allele2": "af",
                     "p.value": "pval_nominal", "BETA": "afc"})
        pos_raw = sus[sus.pip >= POS_PIP].merge(
            nom[["variant_id", "gene", "af", "pval_nominal", "afc"]], on=["variant_id", "gene"])
        pos_raw = qs.parse_colon_ids(pos_raw.dropna(subset=["af"]), args.indel_t)
        cand = pool[pool.orig_gene.isin(pc_genes) & ~pool.variant_id.isin(set(pos_raw.variant_id))]
        pos, m = qs.build_eqtl_context(pos_raw, cand, gene_tpms, annot, args.indel_t,
                                       f"{out}/{ct}", ct)
        for side in ("pos", "neg"):
            qs.set_info_descriptions(f"{out}/{ct}_{side}.vcf", info)
        row = dict(
            cell_type=ct,
            pip90_variants=raw.loc[raw.pip >= POS_PIP, "variantId"].nunique(),
            protein_coding_expressed=sus.loc[sus.pip >= POS_PIP, "variant_id"].nunique(),
            **{vclass: len(pos)}, matched=len(m),
            candidates=len(cand), median_score=float(m.score.median()) if len(m) else np.nan)
        summary.append(row)
        print(f"  {row}", flush=True)

    qs.merge_vcfs(f"{out}/merge.vcf", glob.glob(f"{out}/*_pos.vcf") + glob.glob(f"{out}/*_neg.vcf"))
    qs.merge_matches(f"{out}/matches.tsv", glob.glob(f"{out}/*_matches.tsv"))
    with open(f"{out}/build.json", "w") as fh:
        json.dump({"pip": POS_PIP, "exclude_pip": EXCLUDE_PIP, "indel_t": args.indel_t,
                   "positives": f"TenK10K SuSiE PIP >= pip, protein-coding, measured in {args.expr}",
                   "pool": f"BH q < 0.05 pairs pooled over the {len(cts)} cell types, focal SuSiE genes",
                   "blocklist": "PIP >= exclude_pip in any of 28 cell types",
                   "expression": {"hpa": "HPA rna_immune_cell nTPM",
                                  "onek1k": "OneK1K pseudobulk CPM, Azimuth l2"}[args.expr],
                   "gtf": os.path.basename(args.gtf),
                   "per_cell_type": summary}, fh, indent=2)


if __name__ == "__main__":
    main()
