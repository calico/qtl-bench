# Changelog

## gtex11_v1

First public release. Built from GTEx v11 SuSiE fine-mapping and significant
cis-QTL summary statistics, annotated against GENCODE v48 basic protein-coding
transcripts (readthroughs excluded).

`gtex11_v1.tar.gz` unpacks to `VERSION`, `snp/`, and `indel/`:

| set | pairs | tissues |
|---|---|---|
| `snp/eqtl` | 58,122 | 50 |
| `snp/sqtl` | 58,517 | 50 |
| `snp/paqtl` | 14,823 | 50 |
| `snp/eqtl_gold` | 16,415 | 42 |
| `snp/sqtl_gold` | 19,371 | 47 |
| `snp/paqtl_gold` | 4,548 | 47 |
| `indel/eqtl` | 3,352 * | 50 |
| `indel/sqtl` | 2,857 | 50 |
| `indel/paqtl` | 1,186 | 50 |

\* 3,344 matched pairs: 8 eQTL indel positives exhaust their `(region, Δlen)`
bucket and ship without a negative. Every other set pairs exactly.

SNP sets are built at `--pip 0.9 --exclude_pip 0.01`; indel sets add
`--indel_t 4`. Gold sets use the defaults of `src/gold_set.py`
(`PIP = 1.0`, `NLP ≥ 8`, distance < 500 kb / 100 kb, tissue floor 0.33).

All nine set directories share one `matches.tsv` schema; `mismatch` is identically
zero in the SNP tables, since SNP buckets are exact on alleles.

One known wrinkle: the released SNP positive VCFs were written before positive row
order became canonical, so their equal-PIP ties are ordered arbitrarily. A rebuild
yields the same variants, pairings, and scores, but may permute rows within a tie
group. Key on the ID column rather than row position, and regenerate score files
alongside any rebuilt VCF.

Internal iterations preceding this release were not published.
