#!/usr/bin/env python
"""Convert the TenK10K Phase 1 SuSiE credible sets into the brain_xqtl table schema.

Source: Zenodo 10.5281/zenodo.18221260 (Cuomo et al., "Impact of rare and common genetic
variation on cell type-specific gene expression in human blood"), file susie_summary.zip.
Those per-cell-type files are COMMA-separated despite a .tsv extension, with columns
SNP / Index / PIP / Credible_Set / gene / celltype, where SNP is "chrom:pos:ref:alt"
on GRCh38 without a chr prefix.

Emits {portal}/finemapped_eqtls_pbmc_celltypes_TenK10K_sc.tsv using the same
column names as finemapped_eqtls_brain_celltypes_*.tsv so build_brain_xqtl_sets.py can
consume it. Columns the source does not carry are written as NA:
  rsID, targetGene (symbol), beta, eqtl_beta, Z, pValue, FDR
  non_ref_af  <-- NOTE: the matching score's MAF term needs this; see README_pbmc.md
"""
import glob, os, sys
import numpy as np, pandas as pd

# build_pbmc_finemapped_tsv.py [portal dir]   (default data/tenk10k_pbmc/portal/eqtl)
P = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "../../data/tenk10k_pbmc/portal/eqtl")
SRC = f"{P}/raw/susie_summary"
OUT = f"{P}/finemapped_eqtls_pbmc_celltypes_TenK10K_sc.tsv"

frames = []
for f in sorted(glob.glob(f"{SRC}/*_all_credible_snps.tsv")):
    d = pd.read_csv(f)
    ct = os.path.basename(f).replace("_all_credible_snps.tsv", "")
    p = d.SNP.str.split(":", expand=True)
    d["chrom"] = "chr" + p[0].astype(str)
    d["pos"] = p[1].astype(int)
    d["ref"], d["alt"] = p[2], p[3]
    d["variantId"] = d.chrom + ":" + d.pos.astype(str) + ":" + d.ref + ":" + d.alt
    d["cellType"] = ct
    d["credibleSet"] = ("TenK10K.sc-eQTL." + ct + ".single_context." +
                        d.gene.astype(str) + ".cs95." + d.Credible_Set.astype(str))
    g = d.groupby("credibleSet")
    d["cs95_size"] = g.PIP.transform("size")
    d["is_cs_lead"] = (g.PIP.transform("max") == d.PIP)
    frames.append(d)

a = pd.concat(frames, ignore_index=True)
a = a.rename(columns={"PIP": "pip", "gene": "targetEnsemblID", "celltype": "context"})
a["rsID"] = pd.NA; a["targetGene"] = pd.NA; a["cohort"] = "TenK10K"
for c in ["beta", "eqtl_beta", "Z", "pValue", "FDR", "non_ref_af"]:
    a[c] = np.nan
cols = ["chrom", "pos", "rsID", "variantId", "ref", "alt", "targetGene", "targetEnsemblID",
        "cohort", "context", "cellType", "pip", "beta", "eqtl_beta", "Z", "pValue", "FDR",
        "non_ref_af", "cs95_size", "credibleSet", "is_cs_lead"]
a[cols].to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}\n  {len(a):,} rows, {a.cellType.nunique()} cell types, "
      f"{a.targetEnsemblID.nunique():,} genes, {a.credibleSet.nunique():,} credible sets")
print(f"  PIP>=0.9: {int((a.pip>=0.9).sum()):,}   PIP<0.01: {int((a.pip<0.01).sum()):,}")
