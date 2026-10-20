#!/usr/bin/env python
"""Per-cell-type pseudobulk CPM from OneK1K (Yazar et al. 2022, CELLxGENE collection
dde06e0f-ab3b-46be-96a2-a8082383c4a1), the expression source for PBMC sc-eQTL negative matching.

Sums the raw UMI counts (X, CSR; integer) of every cell of each Azimuth l2 cell type (predicted.celltype.l2,
the annotation TenK10K uses) over all donors, streamed in row
chunks straight from the h5ad, then scales each cell type to counts per million.

  python src/sc/onek1k_pseudobulk.py --h5ad onek1k.h5ad [--out_dir DIR]   ->   pseudobulk_cpm.tsv (Ensembl ID x
                                            cell type) and cell_counts.tsv, the --onek1k_cpm of build_tenk10k_pbmc_eqtl.py
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

CHUNK = 200_000


def cat(f, key):
    """Decode an anndata categorical obs column (codes + categories)."""
    g = f[f"obs/{key}"]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
    return cats[g["codes"][:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="OneK1K h5ad (CELLxGENE collection dde06e0f-ab3b-46be-96a2-a8082383c4a1)")
    ap.add_argument("--out_dir", help="default: the h5ad's directory")
    a = ap.parse_args()
    out = a.out_dir or os.path.dirname(os.path.abspath(a.h5ad))
    with h5py.File(a.h5ad, "r") as f:
        ct = cat(f, "predicted.celltype.l2")  # Azimuth l2, the labels TenK10K uses
        X = f["raw/X"] if "raw" in f else f["X"]
        genes = f["raw/var/_index"][:] if "raw" in f else f["var/_index"][:]
        genes = np.array([g.decode() if isinstance(g, bytes) else g for g in genes])
        codes, labels = pd.factorize(ct)
        indptr = X["indptr"][:]
        n = len(indptr) - 1
        tot = np.zeros((len(labels), len(genes)))
        for s in range(0, n, CHUNK):
            e = min(n, s + CHUNK)
            a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((X["data"][a:b], X["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, len(genes)))
            onehot = sp.csr_matrix((np.ones(e - s), (codes[s:e], np.arange(e - s))), shape=(len(labels), e - s))
            tot += (onehot @ m).toarray()
            print(f"  {e:,}/{n:,} cells", flush=True)
    assert np.allclose(tot, np.round(tot)), "raw/X is not integer counts"
    cpm = pd.DataFrame(tot / tot.sum(1, keepdims=True) * 1e6, index=labels, columns=genes).T
    cpm.index.name = "gene_id"
    cpm.to_csv(f"{out}/pseudobulk_cpm.tsv", sep="\t", float_format="%.4g")
    pd.Series(pd.Series(ct).value_counts(), name="n_cells").to_csv(f"{out}/cell_counts.tsv", sep="\t")
    print(cpm.shape, cpm.columns.tolist())


if __name__ == "__main__":
    main()
