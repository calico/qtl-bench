"""Loaders and plot helpers shared by the GTEx v11 benchmark QC notebooks.

Everything the notebooks need is already in the generated VCFs and match tables,
so the raw parquets are never re-read. Modality-specific INFO fields (`TSSD`/`SD`/
`PD`, `AFC`/`SLOPE`) are normalized on load to `dist` and `effect`, which lets one
notebook loop over all three modalities.
"""

import glob
import os

import numpy as np
import pandas as pd
from scipy import stats as _stats
import matplotlib.pyplot as plt

MODS = ["eqtl", "sqtl", "paqtl"]

DIST_KEY = {"eqtl": "TSSD", "sqtl": "SD", "paqtl": "PD"}                       # VCF INFO
DIST_COL = {"eqtl": "tss_dist", "sqtl": "splice_dist", "paqtl": "pas_dist"}    # matches.tsv
DIST_LABEL = {"eqtl": "TSS distance", "sqtl": "splice site distance",
              "paqtl": "PAS distance"}
EFFECT_KEY = {"eqtl": "AFC", "sqtl": "SLOPE", "paqtl": "SLOPE"}
EFFECT_LABEL = {"eqtl": "|AFC|", "sqtl": "|slope|", "paqtl": "|slope|"}


# --- loading ---------------------------------------------------------------

def parse_vcf(path, info_keys=()):
    """Read a VCF into chrom/pos/variant_id/ref/alt plus the named INFO fields."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            row = {"chrom": c[0], "pos": int(c[1]), "variant_id": c[2],
                   "ref": c[3], "alt": c[4]}
            if info_keys:
                info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
                row.update({k: info.get(k) for k in info_keys})
            rows.append(row)
    return pd.DataFrame(rows)


def read_vcf(path, mod):
    """Read a pos/neg VCF with the modality-specific INFO fields normalized."""
    dist_key, effect_key = DIST_KEY[mod], EFFECT_KEY[mod]
    raw = parse_vcf(path, ["GENE", "TISSUE", "REGION", "PIP", "TPM", "MAF",
                           dist_key, effect_key, "NLP"])
    df = raw[["variant_id", "ref", "alt"]].copy()
    df["gene"], df["tissue"], df["region"] = raw.GENE, raw.TISSUE, raw.REGION
    for name, key in [("pip", "PIP"), ("tpm", "TPM"), ("maf", "MAF"),
                      ("dist", dist_key), ("effect", effect_key), ("nlp", "NLP")]:
        df[name] = pd.to_numeric(raw[key], errors="coerce")
    return df


def load_set(root, mod, labels=("pos", "neg"), tier=""):
    """Concatenate all per-tissue VCFs for one modality, tagged by `label`.

    `tier` selects a filtered sibling set, e.g. tier="_gold" reads {root}/{mod}_gold.
    """
    parts = []
    for label in labels:
        for vcf in sorted(glob.glob(f"{root}/{mod}{tier}/*_{label}.vcf")):
            if os.path.basename(vcf).startswith(("pos_", "neg_")):
                continue  # pos_merge.vcf / neg_merge.vcf
            df = read_vcf(vcf, mod)
            df["label"] = label
            parts.append(df)
    out = pd.concat(parts, ignore_index=True)
    out["ilen"] = out.alt.str.len() - out.ref.str.len()
    return out


def load_matches(root, mod, tier=""):
    """Read the concatenated cross-tissue match table for one modality."""
    return pd.read_csv(f"{root}/{mod}{tier}/matches.tsv", sep="\t")


# --- statistics ------------------------------------------------------------

def mad_zscore(series):
    """Robust z-score: deviation from the median in MAD units."""
    med = series.median()
    mad = (series - med).abs().median()
    return (series - med) / (1.4826 * mad + 1e-12)


def pair_stats(mm, dist_col):
    """Paired correlations and same-gene rate for one matches table.

    `score - mismatch` is the covariate-only score; for SNP tables, whose buckets
    are allele-exact, there is no `mismatch` column and it is just `score`.
    """
    mm = mm.dropna(subset=["neg_variant"])
    out = {"pairs": len(mm)}
    out["r log10(TPM+1)"] = np.corrcoef(np.log10(mm.pos_tpm + 1),
                                        np.log10(mm.neg_tpm + 1))[0, 1]
    # eQTL UTR3/CDS buckets are matched on compartment overlap, not distance, and
    # carry no distance at all.
    d = mm.dropna(subset=[f"pos_{dist_col}", f"neg_{dist_col}"])
    out["r log10(dist+1)"] = np.nan if len(d) < 2 else np.corrcoef(
        np.log10(d[f"pos_{dist_col}"] + 1), np.log10(d[f"neg_{dist_col}"] + 1))[0, 1]
    out["r log10(MAF)"] = np.corrcoef(np.log10(mm.pos_maf.clip(lower=1e-3)),
                                      np.log10(mm.neg_maf.clip(lower=1e-3)))[0, 1]
    out["same gene"] = (mm.pos_gene == mm.neg_gene).mean()
    out["median score"] = (mm.score - mm.get("mismatch", 0)).median()
    return out


def ks_table(df, attrs):
    """Two-sample KS of each tissue against all others, per attribute."""
    rows = []
    for tissue in sorted(df.tissue.unique()):
        g, rest = df[df.tissue == tissue], df[df.tissue != tissue]
        row = {"tissue": tissue}
        for attr in attrs:
            a, b = g[attr].dropna().values, rest[attr].dropna().values
            if len(a) >= 5 and len(b) >= 5:
                ks = _stats.ks_2samp(a, b)
                row[f"ks_{attr}"], row[f"p_{attr}"] = ks.statistic, ks.pvalue
            else:
                row[f"ks_{attr}"], row[f"p_{attr}"] = np.nan, np.nan
        rows.append(row)
    ks_df = pd.DataFrame(rows).set_index("tissue")
    bonf = 0.05 / (len(attrs) * len(ks_df))
    ks_df["max_ks"] = ks_df[[f"ks_{a}" for a in attrs]].max(axis=1)
    ks_df["n_sig"] = (ks_df[[f"p_{a}" for a in attrs]] < bonf).sum(axis=1)
    return ks_df.sort_values("max_ks", ascending=False), bonf


# --- plots -----------------------------------------------------------------

def barh_panels(stats_df, specs, suptitle=None):
    """One horizontal bar chart per (column, title) spec, with a median line."""
    fig, axes = plt.subplots(1, len(specs), figsize=(7 * len(specs), 10))
    for ax, (col, title) in zip(np.atleast_1d(axes), specs):
        s = stats_df[col].sort_values()
        med = s.median()
        s.plot.barh(ax=ax)
        ax.axvline(med, color="k", linestyle="--", lw=1,
                   label=f"median={med:.3f}" if abs(med) < 10 else f"median={med:,.0f}")
        ax.set_title(title)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=8)
    if suptitle:
        fig.suptitle(suptitle, y=1.005, fontweight="bold")
    plt.tight_layout()
    plt.show()


def cdf_overlay(df, col, highlight, xlabel, title, vline=None):
    """CDF of `col` per tissue, with the pooled curve and named tissues emphasized."""
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.sort(df[col].dropna().values)
    ax.plot(x, np.linspace(0, 1, len(x)), color="k", lw=2, label="All tissues", zorder=5)
    cmap = plt.cm.tab20
    for i, tissue in enumerate(sorted(df.tissue.unique())):
        v = np.sort(df.loc[df.tissue == tissue, col].dropna().values)
        ax.plot(v, np.linspace(0, 1, len(v)), color=cmap(i % 20), alpha=0.35, lw=0.8)
    for tissue in highlight:
        v = np.sort(df.loc[df.tissue == tissue, col].dropna().values)
        ax.plot(v, np.linspace(0, 1, len(v)), lw=1.8, label=tissue)
    ax.set_xscale("log")
    if vline is not None:
        ax.axvline(vline, color="gray", linestyle=":", lw=1, label=f"{vline:g}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("CDF")
    ax.set_title(title)
    ax.legend(fontsize=7, loc="lower right")
    plt.tight_layout()
    plt.show()


def mad_heatmap(stats_df, metrics, title):
    """Tissue x metric robust z-scores, oriented so positive = suspicious.

    `metrics` is a list of (column, direction, label); direction "low" flips the
    sign, "both" flags either tail. Returns the z-score frame.
    """
    z_df = pd.DataFrame(index=stats_df.index)
    for col, direction, label in metrics:
        z = mad_zscore(stats_df[col].fillna(stats_df[col].median()))
        z_df[label] = {"low": -z, "both": z.abs()}.get(direction, z)
    z_df["n_flags"] = (z_df > 2).sum(axis=1)
    z_df = z_df.sort_values("n_flags", ascending=False)
    plot_z = z_df.drop(columns="n_flags")

    fig, ax = plt.subplots(figsize=(len(plot_z.columns) * 0.9 + 1, len(plot_z) * 0.27 + 1))
    im = ax.imshow(plot_z.values, aspect="auto", cmap="RdYlGn_r", vmin=-3, vmax=3)
    plt.colorbar(im, ax=ax, label="Robust z-score (+ = suspicious)")
    ax.set_xticks(range(len(plot_z.columns)))
    ax.set_xticklabels(plot_z.columns, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(plot_z)))
    ax.set_yticklabels(plot_z.index, fontsize=7)
    for i in range(len(plot_z)):
        for j in range(len(plot_z.columns)):
            if plot_z.iloc[i, j] > 2:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=False, edgecolor="black", lw=1.5))
    ax.set_title(title)
    plt.tight_layout()
    plt.show()
    return z_df


def overlap_bars(df, a_col, b_col, a_label, b_label, title):
    """Stacked shared / a-only / b-only fractions per tissue, sorted by Jaccard.

    `df` needs columns `tissue`, `shared`, `jaccard`, `a_col` and `b_col`.
    """
    df = df.sort_values("jaccard").reset_index(drop=True)
    union = df.shared + df[a_col] + df[b_col]
    f_shared, f_a = df.shared / union, df[a_col] / union

    fig, ax = plt.subplots(figsize=(16, 6))
    x = np.arange(len(df))
    ax.bar(x, f_shared, color="#2a7", label="shared")
    ax.bar(x, f_a, bottom=f_shared, color="#26c", label=a_label)
    ax.bar(x, df[b_col] / union, bottom=f_shared + f_a, color="#c88", label=b_label)
    for i, jac in enumerate(df.jaccard):
        ax.text(i, 1.01, f"{jac:.2f}", ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(df.tissue, fontsize=8, rotation=45, ha="right")
    ax.set_ylabel("Fraction of union")
    ax.set_ylim(0, 1.12)
    ax.set_title(title)
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


def paired_scatter(ax, x, y, label, title, groups=None):
    """Pos-vs-negative scatter with a y=x line; returns the Pearson r.

    `groups` is an optional list of (mask, color, label) to colour subsets.
    """
    ok = x.notna() & y.notna()
    for sel, color, lbl in (groups or [(ok, "C0", None)]):
        sel = sel & ok
        ax.scatter(x[sel], y[sel], alpha=0.3, s=6, c=color, rasterized=True,
                   label=None if lbl is None else f"{lbl} (n={sel.sum()})")
    lims = [min(x[ok].min(), y[ok].min()), max(x[ok].max(), y[ok].max())]
    ax.plot(lims, lims, "k--", alpha=0.5, lw=1)
    ax.set_xlim(lims); ax.set_ylim(lims); ax.set_aspect("equal")
    ax.set_xlabel(f"positive {label}")
    ax.set_ylabel(f"negative {label}")
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    ax.set_title(f"{title} (r={r:.3f})", fontsize=10)
    if groups:
        ax.legend(fontsize=7, loc="upper left")
    return r
