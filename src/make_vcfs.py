#!/usr/bin/env python
import argparse
import glob
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import pyarrow.parquet as pq
from tqdm import tqdm

"""
make_vcfs.py

Run eqtl_vcfs.py / sqtl_vcfs.py / paqtl_vcfs.py across all GTEx tissues in
parallel and merge per-tissue outputs into study-wide VCFs and matches tables.
"""


################################################################################
# main
################################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Generate positive/negative VCFs across all GTEx tissues.")
    parser.add_argument("--pip", default=0.9, type=float,
                        help="Positive PIP threshold [Default: %(default)s]")
    parser.add_argument("--exclude_pip", default=0.01, type=float,
                        help="Exclude (variant, gene) from negative pool if SuSiE PIP >= this")
    parser.add_argument("--workers", default=4, type=int,
                        help="Parallel tissue workers [Default: %(default)s]")
    parser.add_argument("--qtl", nargs="+", default=["eqtl", "sqtl", "paqtl"],
                        choices=["eqtl", "sqtl", "paqtl"],
                        help="QTL types to process [Default: eqtl sqtl paqtl]")
    parser.add_argument("--susie_dir", default="data/gtex11/portal/GTEx_Analysis_v11_eQTL",
                        help="Directory with eQTL *.SuSiE_summary.parquet files")
    parser.add_argument("--sqtl_dir", default="data/gtex11/portal/GTEx_Analysis_v11_sQTL",
                        help="Directory with sQTL *.SuSiE_summary.parquet files")
    parser.add_argument("--apaqtl_dir", default="data/gtex11/portal/GTEx_Analysis_v11_apaQTL",
                        help="Directory with apaQTL *.SuSiE_summary.parquet files")
    parser.add_argument("--indel_t", default=0, type=int,
                        help="0 = SNPs only; >0 = indels only, up to this length")
    parser.add_argument("--out_root", required=True,
                        help="Root under which {qtltype} set directories are written, "
                             "e.g. data/gtex11/snp")
    parser.add_argument("--cache_dir", default="data/gtex11/portal/cache",
                        help="Shared cache for high_pip_variants.txt / signif_pool.parquet "
                             "[Default: %(default)s]")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    if "eqtl" in args.qtl:
        run_eqtl(args, script_dir)
    if "sqtl" in args.qtl:
        run_sqtl(args, script_dir)
    if "paqtl" in args.qtl:
        run_paqtl(args, script_dir)


################################################################################
# per-QTL-type runners
################################################################################
def prepare_dirs(args, qtl):
    """Make the output set directory and return it with the shared cache paths.

    high_pip_variants.txt and signif_pool.parquet depend only on (qtl, exclude_pip),
    not on pip, out_root, or variant class, so they live in one cache shared across
    versions and across the SNP and indel builds.
    """
    out_dir = f"{args.out_root}/{qtl}"
    cache = f"{args.cache_dir}/{qtl}_excl{int(args.exclude_pip * 10000)}"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(cache, exist_ok=True)
    return out_dir, f"{cache}/high_pip_variants.txt", f"{cache}/signif_pool.parquet"


def write_build_json(out_dir, args, susie_dir):
    """Record the parameters the directory name no longer encodes."""
    with open(f"{out_dir}/build.json", "w") as fh:
        json.dump({"pip": args.pip, "exclude_pip": args.exclude_pip,
                   "indel_t": args.indel_t, "susie_dir": susie_dir}, fh, indent=2)
        fh.write("\n")


def run_eqtl(args, script_dir):
    out_dir, exclude_file, pool_file = prepare_dirs(args, "eqtl")

    susie_files = sorted(glob.glob(f"{args.susie_dir}/*.SuSiE_summary.parquet"))
    if not susie_files:
        raise FileNotFoundError(f"No SuSiE summary parquets in {args.susie_dir}")
    print(f"[eQTL] Found {len(susie_files)} tissue parquets")

    if not os.path.exists(exclude_file):
        collect_high_pip(susie_files, args.exclude_pip, exclude_file, label="eQTL")
    else:
        print(f"[eQTL] Reusing existing {exclude_file}")

    if not os.path.exists(pool_file):
        signif_files = [f.replace("SuSiE_summary", "signif_pairs") for f in susie_files]
        build_signif_pool(signif_files, pool_file, label="eQTL")
    else:
        print(f"[eQTL] Reusing existing {pool_file}")

    eqtl_script = os.path.join(script_dir, "eqtl_vcfs.py")
    jobs = []
    for susie_file in susie_files:
        tissue = os.path.basename(susie_file).split(".v11.")[0]
        out_prefix = f"{out_dir}/{tissue}"
        cmd = (
            f"{sys.executable} {eqtl_script} {susie_file}"
            f" --pos_pip {args.pip:f} --exclude_pip {args.exclude_pip:f}"
            f" --indel_t {args.indel_t}"
            f" --exclude {exclude_file}"
            f" --signif_pool {pool_file}"
            f" -o {out_prefix}"
        )
        jobs.append(cmd)
    exec_par(jobs, args.workers)

    merge_variants(f"{out_dir}/merge.vcf",
                   [f"{out_dir}/*_pos.vcf", f"{out_dir}/*_neg.vcf"])
    # DEPRECATED: split positive/negative merged VCFs, kept only while older
    # two-stage workflows still consume them. Superseded by merge.vcf above;
    # safe to remove once those workflows migrate.
    merge_variants(f"{out_dir}/pos_merge.vcf", f"{out_dir}/*_pos.vcf")
    merge_variants(f"{out_dir}/neg_merge.vcf", f"{out_dir}/*_neg.vcf")
    merge_matches(f"{out_dir}/matches.tsv", f"{out_dir}/*_matches.tsv")

    spec_file = f"{out_dir}/spec_slopes.parquet"
    if not os.path.exists(spec_file):
        build_spec_slopes(susie_files, f"{out_dir}/*_pos.vcf", spec_file)
    else:
        print(f"[eQTL] Reusing existing {spec_file}")

    write_build_json(out_dir, args, args.susie_dir)


def run_sqtl(args, script_dir):
    out_dir, exclude_file, pool_file = prepare_dirs(args, "sqtl")

    susie_files = sorted(glob.glob(f"{args.sqtl_dir}/*.SuSiE_summary.parquet"))
    if not susie_files:
        raise FileNotFoundError(f"No sQTL SuSiE summary parquets in {args.sqtl_dir}")
    print(f"[sQTL] Found {len(susie_files)} tissue parquets")

    if not os.path.exists(exclude_file):
        collect_high_pip(susie_files, args.exclude_pip, exclude_file, label="sQTL")
    else:
        print(f"[sQTL] Reusing existing {exclude_file}")

    if not os.path.exists(pool_file):
        signif_files = [f.replace("SuSiE_summary", "signif_pairs") for f in susie_files]
        build_signif_pool(
            signif_files, pool_file, label="sQTL",
            extra_cols=["group_id"],
        )
    else:
        print(f"[sQTL] Reusing existing {pool_file}")

    sqtl_script = os.path.join(script_dir, "sqtl_vcfs.py")
    jobs = []
    for susie_file in susie_files:
        tissue = os.path.basename(susie_file).split(".v11.")[0]
        out_prefix = f"{out_dir}/{tissue}"
        cmd = (
            f"{sys.executable} {sqtl_script} {susie_file}"
            f" --pos_pip {args.pip:f} --exclude_pip {args.exclude_pip:f}"
            f" --indel_t {args.indel_t}"
            f" --exclude {exclude_file}"
            f" --signif_pool {pool_file}"
            f" -o {out_prefix}"
        )
        jobs.append(cmd)
    exec_par(jobs, args.workers)

    merge_variants(f"{out_dir}/merge.vcf",
                   [f"{out_dir}/*_pos.vcf", f"{out_dir}/*_neg.vcf"])
    # DEPRECATED: split positive/negative merged VCFs, kept only while older
    # two-stage workflows still consume them. Superseded by merge.vcf above;
    # safe to remove once those workflows migrate.
    merge_variants(f"{out_dir}/pos_merge.vcf", f"{out_dir}/*_pos.vcf")
    merge_variants(f"{out_dir}/neg_merge.vcf", f"{out_dir}/*_neg.vcf")
    merge_matches(f"{out_dir}/matches.tsv", f"{out_dir}/*_matches.tsv")

    write_build_json(out_dir, args, args.sqtl_dir)


def run_paqtl(args, script_dir):
    out_dir, exclude_file, pool_file = prepare_dirs(args, "paqtl")

    susie_files = sorted(glob.glob(f"{args.apaqtl_dir}/*.SuSiE_summary.parquet"))
    if not susie_files:
        raise FileNotFoundError(f"No apaQTL SuSiE summary parquets in {args.apaqtl_dir}")
    print(f"[apaQTL] Found {len(susie_files)} tissue parquets")

    if not os.path.exists(exclude_file):
        collect_high_pip(susie_files, args.exclude_pip, exclude_file, label="apaQTL")
    else:
        print(f"[apaQTL] Reusing existing {exclude_file}")

    if not os.path.exists(pool_file):
        signif_files = [f.replace("SuSiE_summary", "signif_pairs") for f in susie_files]
        build_signif_pool(
            signif_files, pool_file, label="apaQTL",
            extra_cols=["group_id"],
        )
    else:
        print(f"[apaQTL] Reusing existing {pool_file}")

    paqtl_script = os.path.join(script_dir, "paqtl_vcfs.py")
    jobs = []
    for susie_file in susie_files:
        tissue = os.path.basename(susie_file).split(".v11.")[0]
        out_prefix = f"{out_dir}/{tissue}"
        cmd = (
            f"{sys.executable} {paqtl_script} {susie_file}"
            f" --pos_pip {args.pip:f} --exclude_pip {args.exclude_pip:f}"
            f" --indel_t {args.indel_t}"
            f" --exclude {exclude_file}"
            f" --signif_pool {pool_file}"
            f" -o {out_prefix}"
        )
        jobs.append(cmd)
    exec_par(jobs, args.workers)

    merge_variants(f"{out_dir}/merge.vcf",
                   [f"{out_dir}/*_pos.vcf", f"{out_dir}/*_neg.vcf"])
    # DEPRECATED: split positive/negative merged VCFs, kept only while older
    # two-stage workflows still consume them. Superseded by merge.vcf above;
    # safe to remove once those workflows migrate.
    merge_variants(f"{out_dir}/pos_merge.vcf", f"{out_dir}/*_pos.vcf")
    merge_variants(f"{out_dir}/neg_merge.vcf", f"{out_dir}/*_neg.vcf")
    merge_matches(f"{out_dir}/matches.tsv", f"{out_dir}/*_matches.tsv")

    write_build_json(out_dir, args, args.apaqtl_dir)


################################################################################
# helpers
################################################################################
def exec_par(cmds, max_proc=None, verbose=True):
    """Run shell commands concurrently, raising if any exits non-zero."""
    def run(cmd):
        if verbose:
            print(cmd, file=sys.stderr)
        return subprocess.run(cmd, shell=True).returncode

    with ThreadPoolExecutor(max_workers=max_proc or len(cmds)) as pool:
        codes = list(pool.map(run, cmds))

    failed = [c for c, code in zip(cmds, codes) if code != 0]
    if failed:
        raise RuntimeError(f"{len(failed)} of {len(cmds)} jobs failed:\n" +
                           "\n".join(failed))


def collect_high_pip(susie_files, pip_threshold, variants_file, label=""):
    """Scan all SuSiE_summary parquets and emit a cross-tissue variant exclusion set.

    variants_file: variant_ids whose max PIP >= threshold in any tissue (one per line).
    """
    high_pip_variants = set()
    desc = f"[{label}] Collecting PIP>={pip_threshold}" if label else f"Collecting PIP>={pip_threshold}"
    for f in tqdm(susie_files, desc=desc):
        df = pq.read_table(
            f, columns=["variant_id", "pip", "biotype"]
        ).to_pandas()
        df = df[(df.biotype == "protein_coding") & (df.pip >= pip_threshold)]
        high_pip_variants.update(df.variant_id.values)

    with open(variants_file, "w") as fh:
        for v in sorted(high_pip_variants):
            fh.write(v + "\n")
    print(f"Wrote {len(high_pip_variants)} high-PIP variants to {variants_file}")


def build_signif_pool(signif_files, pool_file, label="", extra_cols=()):
    """Union all tissues' signif_pairs into a single deduped parquet.

    Each per-tissue *_vcfs.py run draws its negative candidate pool from this file
    instead of the tissue's own signif_pairs, so that low-power tissues are not stuck
    with sparse buckets. af is reduced to its median across the tissues in which the
    pair was significant (af varies slightly by sample subset). Any extra_cols
    (e.g. group_id, start_distance for sQTL/paQTL) are intrinsic to a
    (variant_id, phenotype_id) pair and kept via first() across tissues.
    """
    desc = f"[{label}] Pooling signif_pairs" if label else "Pooling signif_pairs"
    cols = ["variant_id", "phenotype_id", "af", *extra_cols]
    parts = []
    for f in tqdm(signif_files, desc=desc):
        parts.append(pq.read_table(f, columns=cols).to_pandas())
    pooled = pd.concat(parts, ignore_index=True)
    n_raw = len(pooled)
    agg = {"af": "median", **{c: "first" for c in extra_cols}}
    pooled = (
        pooled.groupby(["variant_id", "phenotype_id"], sort=False, as_index=False)
        .agg(agg)
    )
    pooled.to_parquet(pool_file, index=False)
    print(f"Pooled {n_raw:,} -> {len(pooled):,} unique (variant, phenotype) rows in {pool_file}")


def build_spec_slopes(susie_files, pos_vcf_glob, spec_file, label="eQTL"):
    """Assemble a cross-tissue effect-size table for positive (variant, gene) pairs.

    For each variant-gene pair that is a fine-mapped positive in some tissue, gather
    its association slope (with nominal p-value and af) in every tissue where the pair
    is nominally significant, plus the fine-mapped pip/afc where the pair is a SuSiE
    credible-set member. This powers the cross-tissue specificity metrics in
    westminster_eqtl_gtexg.py: slope is dense (~99% of signal cells) and on one scale,
    so it anchors the cross-tissue vectors; the sparser fine-mapped afc is kept as an
    option. Cells absent here are treated as zero-effect downstream.

    Output columns: variant_id, gene (versionless), tissue, slope, pval_nominal, af,
    pip, afc. One row per (variant, gene, tissue) with signal.
    """
    # positive (variant, gene) pairs from the per-tissue positive VCFs
    pos_pairs = set()
    for vcf in sorted(glob.glob(pos_vcf_glob)):
        with open(vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.split("\t")
                info = dict(f.split("=", 1) for f in cols[7].split(";") if "=" in f)
                gene = info.get("GENE", "")
                if gene:
                    pos_pairs.add((cols[2], gene))
    print(f"[{label}] {len(pos_pairs):,} positive (variant, gene) pairs")
    pos_index = pd.MultiIndex.from_tuples(pos_pairs, names=["variant_id", "gene"])

    parts = []
    for susie_file in tqdm(susie_files, desc=f"[{label}] Assembling spec slopes"):
        tissue = os.path.basename(susie_file).split(".v11.")[0]
        signif_file = susie_file.replace("SuSiE_summary", "signif_pairs")

        # association slopes for the positive pairs nominally significant here
        sig = pq.read_table(
            signif_file,
            columns=["phenotype_id", "variant_id", "slope", "pval_nominal", "af"],
        ).to_pandas()
        sig["gene"] = sig.phenotype_id.str.split(".").str[0]
        sig = sig[
            pd.MultiIndex.from_arrays([sig.variant_id, sig.gene]).isin(pos_index)
        ].drop(columns=["phenotype_id"])
        if sig.empty:
            continue

        # fine-mapped pip/afc where the pair is a credible-set member (best PIP)
        sus = pq.read_table(
            susie_file, columns=["variant_id", "phenotype_id", "pip", "afc"]
        ).to_pandas()
        sus["gene"] = sus.phenotype_id.str.split(".").str[0]
        sus = (
            sus.sort_values("pip", ascending=False)
            .drop_duplicates(["variant_id", "gene"])
            .drop(columns=["phenotype_id"])
        )

        sig = sig.merge(sus, on=["variant_id", "gene"], how="left")
        sig.insert(2, "tissue", tissue)
        parts.append(sig)

    spec_df = pd.concat(parts, ignore_index=True)[
        ["variant_id", "gene", "tissue", "slope", "pval_nominal", "af", "pip", "afc"]
    ]
    spec_df.to_parquet(spec_file, index=False)
    print(f"[{label}] Wrote {len(spec_df):,} (pair, tissue) rows to {spec_file}")


def merge_variants(merge_vcf_file, vcf_globs):
    """Merge per-tissue VCFs: write unique variants (by ID) to one sorted VCF.

    vcf_globs may be a single glob string or a list of glob strings.
    """
    if isinstance(vcf_globs, str):
        vcf_globs = [vcf_globs]
    files = sorted(f for g in vcf_globs for f in glob.glob(g))
    if not files:
        print(f"  no files matched {vcf_globs}")
        return

    raw = merge_vcf_file.replace(".vcf", "_raw.vcf")
    with open(raw, "w") as out:
        with open(files[0]) as fh:
            for line in fh:
                if line.startswith("#"):
                    out.write(line)

        seen = set()
        for vcf in files:
            with open(vcf) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    vid = line.split("\t", 3)[2]
                    if vid not in seen:
                        out.write(line)
                        seen.add(vid)

    with open(merge_vcf_file, "w") as out:
        subprocess.run(["bedtools", "sort", "-header", "-i", raw], stdout=out, check=True)
    os.remove(raw)
    print(f"Merged {len(seen)} unique variants into {merge_vcf_file}")


def merge_matches(out_file, matches_glob):
    """Concatenate per-tissue matches tables with a tissue column."""
    dfs = []
    for f in sorted(glob.glob(matches_glob)):
        tissue = os.path.basename(f).split("_matches.tsv")[0]
        df = pd.read_csv(f, sep="\t")
        df.insert(0, "tissue", tissue)
        dfs.append(df)
    if dfs:
        pd.concat(dfs, ignore_index=True).to_csv(out_file, sep="\t", index=False)
        print(f"Merged {len(dfs)} matches tables into {out_file}")


################################################################################
# __main__
################################################################################
if __name__ == "__main__":
    main()
