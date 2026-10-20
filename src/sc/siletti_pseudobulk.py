#!/usr/bin/env python
"""Per-cell-type pseudobulk CPM from the HCA Human Brain Cell Atlas v1.0 (Siletti et al. 2023; CELLxGENE collection
283d65eb-dd53-496d-adb7-7570c7caa443), the expression source for brain sc-eQTL negative matching.

Each benchmark cell type sums the raw UMI counts of the atlas clusters listed for it in siletti_clusters.tsv
(the clusters of the cell-type tracks the benchmark was designed around), over all regions and donors,
streamed in row chunks from the supercluster h5ads (siletti_urls.tsv, CELLxGENE); then CPM.

  python src/sc/siletti_pseudobulk.py --h5ad_dir DIR [--out_dir DIR]   ->   pseudobulk_cpm.tsv (Ensembl ID x
                                             cell type) and cell_counts.tsv, the --siletti_cpm of build_adsp_brain.py
"""
import argparse
import glob
import os

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
CHUNK = 200_000


def cat(f, key):
    g = f[f"obs/{key}"]
    cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in g["categories"][:]])
    return cats[g["codes"][:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad_dir", required=True, help="supercluster h5ads from siletti_urls.tsv")
    ap.add_argument("--out_dir", help="default: --h5ad_dir")
    a = ap.parse_args()
    D, out = a.h5ad_dir, a.out_dir or a.h5ad_dir
    tm = pd.read_csv(f"{HERE}/siletti_clusters.tsv", sep="\t", dtype={"cluster_id": str})
    want = {ct: set(g.cluster_id) for ct, g in tm.groupby("cell_type")}  # cell type -> cluster ids
    owner = {c: ct for ct, cs in want.items() for c in cs}
    cts = sorted(want)
    parts, ncell, seen = [], {ct: 0 for ct in cts}, set()
    for path in sorted(glob.glob(f"{D}/*.h5ad")):
        with h5py.File(path, "r") as f:
            # gene lists differ between supercluster files: sum per file, align on gene ID below
            genes = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var/_index"][:]])
            tot = np.zeros((len(cts), len(genes)))
            cid = cat(f, "cluster_id")
            code = np.array([cts.index(owner[c]) if c in owner else -1 for c in cid])
            seen |= set(cid) & set(owner)
            X = f["X"]
            indptr = X["indptr"][:]
            n = len(indptr) - 1
            for s in range(0, n, CHUNK):
                e = min(n, s + CHUNK)
                keep = np.where(code[s:e] >= 0)[0]
                if not len(keep):
                    continue
                a, b = indptr[s], indptr[e]
                m = sp.csr_matrix((X["data"][a:b], X["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, len(genes)))
                oh = sp.csr_matrix((np.ones(len(keep)), (code[s:e][keep], keep)), shape=(len(cts), e - s))
                tot += (oh @ m).toarray()
            for i, ct in enumerate(cts):
                ncell[ct] += int((code == i).sum())
            parts.append(pd.DataFrame(tot, index=cts, columns=genes))
            print(f"  {path.split('/')[-1]}: {n:,} cells, {(code >= 0).sum():,} in track clusters", flush=True)
    missing = set(owner) - seen
    assert not missing, f"track clusters not found in the downloaded superclusters: {sorted(missing)}"
    tot = pd.concat(parts, axis=1).T.groupby(level=0).sum().T  # a gene absent from a file counts zero there
    assert np.allclose(tot.values, np.round(tot.values)), "X is not integer counts"
    cpm = (tot.div(tot.sum(1), axis=0) * 1e6).T
    cpm.index.name = "gene_id"
    cpm.to_csv(f"{out}/pseudobulk_cpm.tsv", sep="\t", float_format="%.4g")
    pd.Series(ncell, name="n_cells").to_csv(f"{out}/cell_counts.tsv", sep="\t")
    print(pd.Series(ncell))


if __name__ == "__main__":
    main()
