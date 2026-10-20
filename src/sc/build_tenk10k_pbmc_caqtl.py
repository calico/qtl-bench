#!/usr/bin/env python
"""tenk10k_pbmc caqtl: TenK10K PBMC sc-caQTL positive sets (no negatives).

TenK10K releases no per-peak accessibility to match negatives on (the role TPM plays for
eQTLs) and no all-pairs nominal table to draw them from, so the caQTL sets are positives only.
They suit effect-size evaluation (e.g. Spearman of predicted vs SLOPE), not classification.

  positives   SuSiE PIP >= 0.9 (all_celltypes_caQTL_susie_annotated/{ct}_susie_annotated.csv),
              SNPs or indels <= 4 bp (--indel_t), one row per variant (its highest-PIP peak);
              PKD = distance to that peak (0 inside), SLOPE and NLP from the SuSiE CSV

  python src/sc/build_tenk10k_pbmc_caqtl.py --indel_t 0    # -> data/tenk10k_pbmc/snp/caqtl/
  python src/sc/build_tenk10k_pbmc_caqtl.py --indel_t 4    # -> data/tenk10k_pbmc/indel/caqtl/
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

import qb_sets as qs

PORTAL = f"{qs.REPO}/data/tenk10k_pbmc/portal/caqtl"  # raw/ from fetch_tenk10k_pbmc.sh
POS_PIP = 0.9


def parse_peaks(peaks):
    p = pd.Series(peaks).str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
    return p[0].values, p[1].astype(np.int64).values, p[2].astype(np.int64).values


def peak_dist(pos, start, end):
    """0 inside the peak, else distance to its nearest edge."""
    return np.where(pos < start, start - pos, np.where(pos > end, pos - end, 0)).astype(float)


def load_susie(susie, ct):
    d = pd.read_csv(f"{susie}/{ct}_susie_annotated.csv")
    d["variant_id"] = "chr" + d.SNP
    return d.rename(columns={"PIP": "pip"})


def write_vcf(path, df, tissue):
    header = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=PEAK,Number=1,Type=String,Description="caQTL peak (the variant\'s highest-PIP peak)">',
        '##INFO=<ID=TISSUE,Number=1,Type=String,Description="TenK10K cell type (Azimuth l2)">',
        '##INFO=<ID=PIP,Number=1,Type=Float,Description="SuSiE posterior inclusion probability">',
        '##INFO=<ID=PKD,Number=1,Type=Integer,Description="Distance to PEAK (0 inside)">',
        '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Alt allele frequency">',
        '##INFO=<ID=NLP,Number=1,Type=Float,Description="Negative log10 nominal p-value">',
        '##INFO=<ID=SLOPE,Number=1,Type=Float,Description="tensorQTL slope (alt allele)">',
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]),
    ]
    with open(path, "w") as out:
        out.write("\n".join(header) + "\n")
        for v in df.itertuples():
            nlp = "." if pd.isna(v.pval_nominal) else f"{min(-np.log10(max(v.pval_nominal, 1e-300)), 300):.4f}"
            slope = "." if pd.isna(v.slope) else f"{v.slope:.4f}"
            info = (f"PEAK={v.peak};TISSUE={tissue};PIP={v.pip:.4f};PKD={int(v.pk_dist)};"
                    f"MAF={v.maf:.4f};AF={v.af:.4f};NLP={nlp};SLOPE={slope}")
            out.write("\t".join([v.chrom, str(v.position), v.variant_id, v.ref, v.alt,
                                 ".", ".", info]) + "\n")


def main():
    qs.pin_hash_seed()
    ap = argparse.ArgumentParser()
    ap.add_argument("--indel_t", type=int, default=0)
    ap.add_argument("--portal", default=PORTAL, help="TenK10K caQTL inputs (raw/)")
    ap.add_argument("--out", help="output set directory (default data/tenk10k_pbmc/{snp,indel}/caqtl)")
    args = ap.parse_args()
    vclass = "snp" if args.indel_t == 0 else "indel"
    susie_dir = f"{args.portal}/raw/all_celltypes_caQTL_susie_annotated"
    out = args.out or f"{qs.REPO}/data/tenk10k_pbmc/{vclass}/caqtl"
    os.makedirs(out, exist_ok=True)
    cts = sorted(os.path.basename(f)[:-len("_susie_annotated.csv")]
                 for f in glob.glob(f"{susie_dir}/*_susie_annotated.csv"))

    summary = []
    for ct in cts:
        raw = load_susie(susie_dir, ct)
        pos = qs.parse_colon_ids(raw[raw.pip >= POS_PIP], args.indel_t)
        pos = pos.sort_values(["pip", "variant_id"], ascending=[False, True], kind="stable")
        pos = pos.drop_duplicates("variant_id")
        _, s, e = parse_peaks(pos.peak.values)
        pos["pk_dist"] = peak_dist(pos.position.values, s, e)
        write_vcf(f"{out}/{ct}_pos.vcf", pos, ct)
        row = dict(cell_type=ct, pip90_variants=raw.loc[raw.pip >= POS_PIP, "variant_id"].nunique(),
                   **{vclass: len(pos)}, in_peak=int((pos.pk_dist == 0).sum()))
        summary.append(row)
        print(f"  {row}", flush=True)

    qs.merge_vcfs(f"{out}/merge.vcf", glob.glob(f"{out}/*_pos.vcf"))
    with open(f"{out}/build.json", "w") as fh:
        json.dump({"pip": POS_PIP, "indel_t": args.indel_t,
                   "positives": "SuSiE PIP >= pip, highest-PIP peak per variant",
                   "negatives": "none (TenK10K releases no per-peak accessibility or all-pairs table)",
                   "per_cell_type": summary}, fh, indent=2)


if __name__ == "__main__":
    main()
