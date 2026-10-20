#!/usr/bin/env python
"""One pass over TenK10K common_all_cis_pvalues_100kb.zip per cell type.

Writes to {portal}/nominal/:
  {ct}_pos_beta.tsv     fine-mapped positives (PIP >= 0.9) with SAIGE-QTL BETA, SE, p, AF
                        -- the effect sizes finemapped_eqtls_pbmc_celltypes_TenK10K_sc.tsv
                        leaves as NA
  {ct}_bh_pool.tsv.gz   negative pool: variant-gene pairs passing per-gene BH (q < --fdr),
                        BH applied within each gene's cis tests (decided 2026-09-23)

Allele2 is the alt allele (MarkerID is CHR:POS:Allele1:Allele2), so BETA and AF_Allele2 are
alt-allele quantities, matching the positives' variantId orientation. The cis window is gene
body +/- 100 kb, not TSS-centred.
"""
import argparse
import os
import zipfile

import numpy as np
import pandas as pd

COLS = ["MarkerID", "AF_Allele2", "BETA", "SE", "p.value", "gene", "is_snp"]

ap = argparse.ArgumentParser()
ap.add_argument("--celltypes", nargs="+",
                default=["CD4_Naive", "NK", "CD4_TCM", "CD8_TEM", "B_naive", "CD14_Mono", "CD16_Mono"])
ap.add_argument("--fdr", type=float, default=0.05)
ap.add_argument("--pip", type=float, default=0.9)
ap.add_argument("--portal", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "../../data/tenk10k_pbmc/portal/eqtl"))
a = ap.parse_args()
E = a.portal

fm = pd.read_csv(f"{E}/finemapped_eqtls_pbmc_celltypes_TenK10K_sc.tsv", sep="\t",
                 usecols=["variantId", "targetEnsemblID", "cellType", "pip"])
z = zipfile.ZipFile(f"{E}/common_all_cis_pvalues_100kb.zip")


def bh(p):
    """Benjamini-Hochberg q-values for one gene's tests."""
    p = np.asarray(p)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


for ct in a.celltypes:
    d = pd.read_csv(z.open(f"{ct}_common_all_cis_raw_pvalues.tsv"), sep="\t", usecols=COLS,
                    dtype={"MarkerID": str, "gene": str})
    d["variant"] = "chr" + d.MarkerID
    d["q"] = d.groupby("gene")["p.value"].transform(bh)
    keep = ["variant", "gene", "BETA", "SE", "p.value", "q", "AF_Allele2", "is_snp"]

    pos = fm[(fm.cellType == ct) & (fm.pip >= a.pip)]
    pos = pos.merge(d[keep], left_on=["variantId", "targetEnsemblID"],
                    right_on=["variant", "gene"], how="left")
    pos.drop(columns=["variant", "gene"]).to_csv(f"{E}/nominal/{ct}_pos_beta.tsv", sep="\t", index=False)

    pool = d.loc[d.q < a.fdr, keep]
    pool.to_csv(f"{E}/nominal/{ct}_bh_pool.tsv.gz", sep="\t", index=False)
    print(f"{ct}: {len(d):,} tests, {d.gene.nunique():,} genes | positives {len(pos):,} "
          f"(beta found {pos.BETA.notna().mean():.1%}) | BH pool {len(pool):,} pairs, "
          f"{pool.gene.nunique():,} genes", flush=True)
