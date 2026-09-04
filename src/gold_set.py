#!/usr/bin/env python
import argparse
import glob
import json
import os

import pandas as pd

import make_vcfs
import qc_lib

"""
gold_set.py

Derive a high-confidence subset of a benchmark set directory. Positives are
filtered on fine-mapping confidence (PIP), nominal significance (NLP), and
distance to the cognate feature; each rejected positive's matched negative is
dropped with it so the pairing stays one-to-one; tissues left with too few
positives are dropped whole. The output is the same file inventory as the input,
so every consumer works unchanged. See gold_criteria.md for the rationale.
"""

DIST_CAP = {"eqtl": 500_000, "sqtl": 100_000, "paqtl": 100_000}


################################################################################
# main
################################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Filter a QTL benchmark set down to high-confidence positives.")
    parser.add_argument("set_dir",
                        help="Source set directory, e.g. data/gtex11/snp/eqtl")
    parser.add_argument("-o", dest="out_dir", default=None,
                        help="Output directory [Default: {set_dir}_gold]")
    parser.add_argument("--qtl", default=None, choices=list(DIST_CAP),
                        help="Modality [Default: basename of set_dir]")
    parser.add_argument("--pip", default=1.0, type=float,
                        help="Minimum SuSiE PIP [Default: %(default)s]")
    parser.add_argument("--nlp", default=8.0, type=float,
                        help="Minimum -log10 nominal p; positives without one are "
                             "always dropped [Default: %(default)s]")
    parser.add_argument("--dist", default=None, type=int,
                        help=f"Maximum feature distance [Default: {DIST_CAP}]")
    parser.add_argument("--tissue_frac", default=0.33, type=float,
                        help="Drop tissues below this fraction of the median surviving "
                             "positive count [Default: %(default)s]")
    args = parser.parse_args()

    qtl = args.qtl or os.path.basename(args.set_dir.rstrip("/"))
    if qtl not in DIST_CAP:
        raise ValueError(f"Cannot infer modality from {args.set_dir}; pass --qtl")
    dist_key = qc_lib.DIST_KEY[qtl]
    dist_cap = args.dist or DIST_CAP[qtl]
    out_dir = args.out_dir or f"{args.set_dir.rstrip('/')}_gold"
    os.makedirs(out_dir, exist_ok=True)

    # filter every tissue, then drop the tissues that came out too small
    kept = {t: filter_tissue(args.set_dir, t, dist_key, dist_cap, args)
            for t in tissues(args.set_dir)}
    counts = pd.Series({t: len(v[0][1]) for t, v in kept.items()})
    floor = args.tissue_frac * counts.median()
    dropped = sorted(counts.index[counts < floor])
    print(f"[{qtl}] median {counts.median():.0f} positives/tissue, floor {floor:.0f}; "
          f"dropping {len(dropped)} tissues: {', '.join(dropped)}")

    for tissue in sorted(set(kept) - set(dropped)):
        write_tissue(out_dir, tissue, *kept[tissue])

    make_vcfs.merge_variants(f"{out_dir}/merge.vcf",
                             [f"{out_dir}/*_pos.vcf", f"{out_dir}/*_neg.vcf"])
    make_vcfs.merge_matches(f"{out_dir}/matches.tsv", f"{out_dir}/*_matches.tsv")

    n_pos = int(counts.drop(dropped).sum())
    with open(f"{out_dir}/gold.json", "w") as fh:
        json.dump({"source": args.set_dir, "qtl": qtl,
                   "pip": args.pip, "nlp": args.nlp, "dist": dist_cap,
                   "tissue_frac": args.tissue_frac,
                   "tissues": len(counts) - len(dropped), "dropped_tissues": dropped,
                   "positives": n_pos}, fh, indent=2)
        fh.write("\n")
    print(f"[{qtl}] {n_pos:,} positives across {len(counts) - len(dropped)} tissues "
          f"in {out_dir}")


################################################################################
# helpers
################################################################################
def tissues(set_dir):
    return sorted(os.path.basename(f)[: -len("_pos.vcf")]
                  for f in glob.glob(f"{set_dir}/*_pos.vcf")
                  if not os.path.basename(f).startswith(("pos_", "neg_")))


def read_vcf_lines(path):
    """Split a VCF into its header lines and (variant_id, INFO dict, line) records."""
    header, records = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            c = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
            records.append((c[2], info, line))
    return header, records


def keep_positive(info, dist_key, dist_cap, args):
    """Confidence, significance, and range. A missing distance passes: eQTL exonic
    positives carry no TSSD and are inside the gene by construction."""
    if float(info.get("PIP", 0)) < args.pip:
        return False
    nlp = info.get("NLP", ".")
    if nlp == "." or float(nlp) < args.nlp:
        return False
    dist = info.get(dist_key, ".")
    return dist == "." or abs(float(dist)) < dist_cap


def filter_tissue(set_dir, tissue, dist_key, dist_cap, args):
    """Return (pos, neg, matches) restricted to the surviving positives and their
    matched negatives. Unmatched positives keep their empty matches row."""
    pos_header, pos_recs = read_vcf_lines(f"{set_dir}/{tissue}_pos.vcf")
    neg_header, neg_recs = read_vcf_lines(f"{set_dir}/{tissue}_neg.vcf")
    matches = pd.read_csv(f"{set_dir}/{tissue}_matches.tsv", sep="\t")

    pos_ids = {vid for vid, info, _ in pos_recs
               if keep_positive(info, dist_key, dist_cap, args)}
    matches = matches[matches.pos_variant.isin(pos_ids)]
    neg_ids = set(matches.neg_variant.dropna())

    pos = (pos_header, [ln for vid, _, ln in pos_recs if vid in pos_ids])
    neg = (neg_header, [ln for vid, _, ln in neg_recs if vid in neg_ids])
    return pos, neg, matches


def write_tissue(out_dir, tissue, pos, neg, matches):
    for label, (header, lines) in [("pos", pos), ("neg", neg)]:
        with open(f"{out_dir}/{tissue}_{label}.vcf", "w") as fh:
            fh.writelines(header)
            fh.writelines(lines)
    matches.to_csv(f"{out_dir}/{tissue}_matches.tsv", sep="\t", index=False)


################################################################################
# __main__
################################################################################
if __name__ == "__main__":
    main()
