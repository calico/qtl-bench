# qtl-bench

Benchmark sets for evaluating sequence-to-function models on regulatory variant
effect prediction, built from statistically fine-mapped QTLs and their matched
negatives.

Three studies are built:

| study | release | QTLs | source |
|---|---|---|---|
| `gtex11` | `gtex11_v1` | bulk eQTL, sQTL, apaQTL · 50 tissues | GTEx v11 |
| `adsp_brain` | `adsp_brain_v1` | single-nucleus eQTL · 6 brain cell types, 3 ROSMAP cohorts | ADSP FunGen xQTL (ROSMAP) |
| `tenk10k_pbmc` | `tenk10k_pbmc_v1` | single-cell eQTL (28 cell types) and caQTL positives (26) · PBMC | TenK10K Phase 1 |

The single-cell studies are described in [Single-cell studies](#single-cell-studies);
everything else below is GTEx. **GTEx v11** (`gtex11`) covers three molecular
phenotypes across all 50 tissues:

- **eQTL** — cis gene expression
- **sQTL** — cis splicing (LeafCutter intron-cluster usage)
- **apaQTL** — cis alternative polyadenylation (3′ UTR usage ratios)

For each modality and tissue the pipeline emits a positive set (SuSiE PIP ≥ 0.9) and
a one-to-one matched negative set drawn from the same modality's pool of
nominally-significant cis variants, controlled for expression of the linked gene,
distance to the relevant regulatory feature, minor allele frequency, and alleles.
The matching is the point: a model that only learns "causal variants sit near
promoters" cannot beat a negative set that sits equally near promoters.

Separate parallel sets are built for SNPs and for indels, and a **gold** subset of
each restricts to positives whose causal identity and measured effect are both
unambiguous.

## Get the data

```
./fetch.sh gtex11_v1
```

which unpacks into `data/gtex11/`:

```
data/gtex11/
  snp/{eqtl,sqtl,paqtl}/          full SNP sets, PIP >= 0.9
  snp/{eqtl,sqtl,paqtl}_gold/     high-confidence subsets
  indel/{eqtl,sqtl,paqtl}/        indel sets, <= 4 bp
  VERSION
```

Every one of those nine directories has the same file inventory, so any consumer
that takes one path works with all of them.

| File | Contents |
|---|---|
| `{tissue}_pos.vcf` | Fine-mapped positive variants for that tissue, sorted by PIP descending |
| `{tissue}_neg.vcf` | One score-matched negative per positive, where a match exists |
| `{tissue}_matches.tsv` | Paired rows with both members' TPM / MAF / distance / PIP / score |
| `merge.vcf` | Positives and negatives together, deduplicated and sorted — what model scoring reads |
| `matches.tsv` | All per-tissue match tables concatenated, with a `tissue` column |
| `build.json` | `pip`, `exclude_pip`, `indel_t`, `susie_dir` for this build |
| `gold.json` | Gold directories only: the criteria and resulting counts |
| `spec_slopes.parquet` | eQTL only: cross-tissue slope/PIP/AFC table for positive (variant, gene) pairs |
| `pos_merge.vcf` / `neg_merge.vcf` | Deprecated split merges, superseded by `merge.vcf`; not emitted for gold sets |

Set sizes at PIP ≥ 0.9 (positive–negative pairs, summed over tissues):

| | full | gold |
|---|---|---|
| eQTL | 58,122 · 50 tissues | 16,415 · 42 |
| sQTL | 58,517 · 50 | 19,371 · 47 |
| apaQTL | 14,823 · 50 | 4,548 · 47 |

Cross-tissue unique positive variants: ~23.8k eQTL, ~20.7k sQTL, ~5.0k apaQTL.

### VCF INFO fields

Downstream evaluation reads these directly — no rejoin to the source parquets.

**eQTL:**

| Field | Description |
|---|---|
| `GENE` | Ensembl gene ID (versionless) for the eQTL target |
| `TISSUE` | GTEx tissue name |
| `PIP` | SuSiE posterior inclusion probability (0 on negatives) |
| `REGION` | `UTR3`, `CDS`, or `TSS` |
| `TSSD` | Distance to nearest TSS of the eQTL gene (`.` for UTR3/CDS) |
| `TPM` | GTEx median TPM of the gene in the tissue |
| `MAF` | Minor allele frequency |
| `AFC` | SuSiE allelic fold change (log2, alt vs ref); positives only |
| `NLP` | −log10 nominal p-value from tensorQTL (capped at 300); positives only |

**sQTL and apaQTL:**

| Field | Description |
|---|---|
| `GENE` | SuSiE-associated gene for positives; genome-wide nearest splice-site / PAS gene for negatives |
| `TISSUE` | GTEx tissue |
| `PIP` | SuSiE PIP (0 on negatives) |
| `SD` (sQTL) / `PD` (apaQTL) | Distance to the relevant feature |
| `TPM`, `MAF` | as above |
| `NLP` | −log10 nominal p-value; positives only |
| `SLOPE` | tensorQTL per-allele slope; positives only |

## Gold sets

`src/gold_set.py` derives a high-confidence subset of any set directory, for faster
evaluation on cleaner positives. It filters on three INFO fields already present,
drops each rejected positive's matched negative so the pairing stays one-to-one, then
drops tissues left with too few positives. Output is the same inventory under a
`_gold` sibling, plus `gold.json`.

| criterion | eQTL | sQTL | apaQTL |
|---|---|---|---|
| distance cap | `TSSD` < 500 kb | `SD` < 100 kb | `PD` < 100 kb |
| significance | `NLP` present and ≥ 8 | same | same |
| confidence | `PIP` = 1.0 | same | same |
| tissue floor | ≥ 0.33 × the modality's median surviving count | same | same |

eQTL exonic positives (`REGION` in `CDS`/`UTR3`) carry no `TSSD` and are inside the
gene by construction, so they pass the distance cap unconditionally.

The net effect is a 3–3.5× reduction in set size. See
[docs/gold_criteria.md](docs/gold_criteria.md) for the evidence behind each
threshold, including why there is deliberately no effect-size floor.

## Indel sets

`indel/` holds insertion/deletion counterparts built by the same pipeline with
`--indel_t 4`. Construction is identical apart from the variant-class filter and
[allele matching](#allele-matching). Every record is a pure left-anchored indel of
≤ 4 bp, every matched pair shares `ilen` exactly, and `matches.tsv` carries a
`mismatch` column; `score − mismatch` is the covariate-only score, directly
comparable to the SNP sets.

**Why 4 bp.** Longer indels displace the model's output bins more, complicating
scoring, and raising the cap to 8 bp adds only ~10% more positives. The `ilen` bucket
must also stay populated, and 5–8 bp events are an order of magnitude rarer than
1–4 bp ones.

| | unique positives | tissue-level positives | matched | exact alleles |
|---|---|---|---|---|
| `eqtl` | 1,392 | 3,352 | 99.8% | 70.1% |
| `sqtl` | 1,105 | 2,857 | 100.0% | 69.8% |
| `paqtl` | 349 | 1,186 | 100.0% | 54.0% |

eQTL match rate by region: TSS 100.0% (2,768), UTR3 100.0% (363), CDS 96.4% (221).

The allele relaxation improves covariate matching as well as match rate: positives
that would otherwise take a poor exact-allele negative take a better-controlled one
for a penalty of one or two bases.

| median absolute difference | eQTL | sQTL | apaQTL |
|---|---|---|---|
| Δlog10 TPM | 0.138 → **0.116** | 0.151 → **0.120** | 0.221 → **0.155** |
| Δlog2 MAF | 0.292 → **0.255** | 0.384 → **0.322** | 0.535 → **0.441** |
| Δlog10 distance | 0.116 → **0.085** | 0.146 → **0.109** | 0.452 → **0.205** |

(left of each arrow: exact-allele bucketing; right: `ilen` bucketing with the
mismatch penalty). Matching is still looser than for SNPs, which reach
0.016 / 0.033 / 0.011 on the same three measures for eQTL — indel buckets simply hold
fewer candidates.

**Scoring.** Length-changing variants need reference predictions spliced at the
variant bin rather than duplicated across the left/right compensation shifts;
in `hound_snp.py` this is `--indel_stitch`.

## Single-cell studies

`adsp_brain` and `tenk10k_pbmc` apply the GTEx recipe to cell-type-resolved QTLs. A "tissue"
is a cell type. eQTL sets are built exactly as GTEx eQTLs are: SuSiE PIP ≥ 0.9 positives, and
one-to-one negatives drawn from a pool of nominally significant cis variants, matched on
expression, TSS distance, MAF and alleles with the same score and greedy assignment
(`src/eqtl_vcfs.py` functions, called from `src/sc/`). caQTL sets are positives only.

```
./fetch.sh adsp_brain_v1      # -> data/adsp_brain/
./fetch.sh tenk10k_pbmc_v1    # -> data/tenk10k_pbmc/

data/adsp_brain/
  snp/{eqtl,eqtl_rosmap_mit,eqtl_rosmap_cuimc1}/  indel/{...}/  VERSION
data/tenk10k_pbmc/
  snp/{eqtl,caqtl}/  indel/{eqtl,caqtl}/  VERSION
```

eQTL set directories have the GTEx file inventory (`{cell_type}_pos.vcf`, `_neg.vcf`,
`_matches.tsv`, `merge.vcf`, `matches.tsv`, `build.json`). caQTL set directories hold only
`{cell_type}_pos.vcf`, `merge.vcf` (the positives) and `build.json`. There are no `_gold`
sets, `spec_slopes.parquet` or deprecated `pos_merge`/`neg_merge` files.

| set | pairs (eQTL) / positives (caQTL) | cell types | unique positives |
|---|---|---|---|
| `adsp_brain/snp/eqtl` | 3,842 | 6 | 3,115 |
| `adsp_brain/snp/eqtl_rosmap_mit` | 2,850 | 6 | 2,386 |
| `adsp_brain/snp/eqtl_rosmap_cuimc1` | 2,284 | 6 | 1,856 |
| `adsp_brain/indel/eqtl` | 692 | 6 | 639 |
| `adsp_brain/indel/eqtl_rosmap_mit` | 665 | 6 | 621 |
| `adsp_brain/indel/eqtl_rosmap_cuimc1` | 446 | 6 | 409 |
| `tenk10k_pbmc/snp/eqtl` | 13,544 | 28 | 6,350 |
| `tenk10k_pbmc/indel/eqtl` | 2,969 | 28 | 1,791 |
| `tenk10k_pbmc/snp/caqtl` | 38,672 positives | 26 | 25,493 |
| `tenk10k_pbmc/indel/caqtl` | 22,200 positives | 26 | 18,577 |

Every eQTL positive has a matched negative.

Sources and cell types:
- **adsp_brain:** ADSP Functional Genomics xQTL Atlas (NIAGADS NG00184.v1), ROSMAP snRNA
  eQTLs. `eqtl` is the combined CUIMC1+CUIMC2+MIT cohort; `eqtl_rosmap_mit` and
  `eqtl_rosmap_cuimc1` are the MIT and CUIMC1 cohorts fine-mapped alone. The cohorts share
  donors, so they are not independent replicates. Fine-mapping on different donor sets still
  picks different variants: only 41% (MIT) and 58% (CUIMC1) of their positives are positives in
  the combined set. Cell types: `excitatory_neuron`, `inhibitory_neuron`, `oligodendrocyte`,
  `astrocyte`, `OPC`, `microglia`.
- **tenk10k_pbmc eQTL:** TenK10K Phase 1 scRNA (Zenodo 10.5281/zenodo.18221260), all 28
  Azimuth l2 cell types, from `CD4_TCM` (1,525 SNP pairs) to `ILC` (21).
- **tenk10k_pbmc caQTL:** TenK10K Phase 1 multiome (Hugging Face `anglixue/TenK10K_multiome`),
  all 26 released cell types, from `CD14_Mono` (6,213 SNP positives) to `ASDC` (136).

### Differences from gtex11

- **Variant IDs** are `chrN:pos:ref:alt` (GRCh38), as the sources release them, not
  GTEx's `chrN_pos_ref_alt_b38`. Key on the ID column.
- **eQTL INFO** has the GTEx keys, with these meanings (the VCF headers say the same):

  | Field | `adsp_brain` | `tenk10k_pbmc` |
  |---|---|---|
  | `TISSUE` | FunGen snuc-eQTL cell type | TenK10K Azimuth l2 cell type |
  | `TPM` | HCA Human Brain Cell Atlas v1.0 (Siletti et al. 2023) pseudobulk CPM over the cell type's atlas clusters | OneK1K (Yazar et al. 2022) pseudobulk CPM, Azimuth l2 |
  | `AFC` | FunGen SuSiE `conditional_effect` (alt vs ref) | SAIGE-QTL `BETA` (alt allele) |
  | `NLP` | −log10 FunGen `qtl_p_val` | −log10 SAIGE-QTL nominal p |

  `AFC` is therefore not a log2 allelic fold change on these sets; compare it within a
  study only.
- **caQTL sets have no negatives.** Matching needs per-peak accessibility (the role TPM plays
  for eQTLs) and a pool of nominally significant non-causal variants. TenK10K releases neither,
  so the caQTL sets support effect-size evaluation (e.g. Spearman of predicted against `SLOPE`)
  but not classification. caQTL INFO: `PEAK` (the variant's highest-PIP peak), `TISSUE`, `PIP`,
  `PKD` (distance to `PEAK`, 0 inside), `MAF`, `AF`, `NLP` and `SLOPE` (tensorQTL, alt allele).
- **No gold sets yet.** `gold_set.py` is GTEx-calibrated, and the single-cell sets are
  small enough to evaluate whole.

### Construction

| | `adsp_brain` eQTL | `tenk10k_pbmc` eQTL | `tenk10k_pbmc` caQTL |
|---|---|---|---|
| positives | FunGen `single_context_finemapping_all` PIP ≥ 0.9, protein-coding gene in the HCA atlas | SuSiE PIP ≥ 0.9, protein-coding gene in OneK1K; AF, p and BETA from the 100 kb SAIGE-QTL tests | SuSiE PIP ≥ 0.9, highest-PIP peak per variant |
| pool | the cohort's `bh_significant` pairs of the six cell types (AF = median), restricted to the focal cell type's fine-mapped genes | per-gene BH q < 0.05 pairs of all 28 cell types, focal fine-mapped genes | — |
| blocklist | PIP ≥ 0.01 in any cell type of any ROSMAP cohort | PIP ≥ 0.01 in any of the 28 cell types | — |
| negative feature | nearest TSS (GTEx rule) | nearest TSS | — |
| buckets | (region, ref, alt) / (region, Δlen) | same | — |

Rebuild:

```
src/sc/fetch_adsp_brain.sh                      # NIAGADS NG00184 -> data/adsp_brain/portal/raw/
src/sc/fetch_tenk10k_pbmc.sh all                # Zenodo + HF -> data/tenk10k_pbmc/portal/{eqtl,caqtl}/
python src/sc/build_pbmc_finemapped_tsv.py      # TenK10K susie_summary -> fine-mapping table
python src/sc/extract_pbmc_eqtl_nominal.py --celltypes <all 28>   # positives' BETA/AF/p and the BH pool
python src/sc/siletti_pseudobulk.py --h5ad_dir SIL    # HCA atlas h5ads (siletti_urls.tsv) -> pseudobulk_cpm.tsv
python src/sc/onek1k_pseudobulk.py --h5ad onek1k.h5ad  # CELLxGENE dde06e0f-... -> pseudobulk_cpm.tsv
python src/sc/build_adsp_brain.py --indel_t {0,4} --cohort {ROSMAP_CUIMC1_2_MIT,ROSMAP_MIT,ROSMAP_CUIMC1} \
    --siletti_cpm SIL/pseudobulk_cpm.tsv
python src/sc/build_tenk10k_pbmc_eqtl.py --indel_t {0,4} --onek1k_cpm ONEK1K/pseudobulk_cpm.tsv
python src/sc/build_tenk10k_pbmc_caqtl.py --indel_t {0,4}
```

Each brain cell type's pseudobulk sums the atlas clusters listed in `src/sc/siletti_clusters.tsv`
(153 clusters). The builders read `data/gtex11/portal/gencode48_basic_nort_protein.gtf`
(`portal.mk`), and every input is a public download. The builders pin `PYTHONHASHSEED=0` (set
order otherwise breaks a few near-ties differently run to run), so rebuilds reproduce the
released files byte for byte.

## Rebuilding from source

Everything below is only needed to regenerate the sets; to use them, `fetch.sh` is
enough.

### Inputs

```
make -f portal.mk all        # ~13 GB, download-bound
```

downloads into `data/gtex11/portal/`:

| Source | Purpose |
|---|---|
| GTEx v11 SuSiE fine-mapping tarballs | Per-tissue variant–phenotype PIPs, allelic fold change (eQTL), credible-set IDs, biotype |
| GTEx v11 significant cis-QTL tarballs | Per-tissue `signif_pairs` — the wider pool of nominally significant (variant, phenotype) pairs; source of negative candidates and per-positive p-values/slopes |
| GTEx v11 median TPM GCT | Per-(tissue, gene) expression for matching and gene-existence filtering |
| GENCODE v48 basic annotation | CDS / UTR intervals, TSS positions, MANE_Select / GENCODE_Primary tags, splice sites, PAS sites |

The annotation target filters GENCODE v48 basic down to protein-coding transcripts
and strips readthroughs, which would otherwise contaminate UTR3 annotation with
spuriously extended 3′ UTRs and pollute the splice-site and PAS sets.

### Build

```
# SNPs, all modalities, all tissues
python src/make_vcfs.py --pip 0.9 --workers 12 --out_root data/gtex11/snp

# indels
python src/make_vcfs.py --pip 0.9 --workers 4 --indel_t 4 --out_root data/gtex11/indel

# gold subsets
for q in eqtl sqtl paqtl; do python src/gold_set.py data/gtex11/snp/$q; done
```

A few minutes per modality at 12 workers; the heaviest single step is GTF parsing in
each per-tissue process. `gold_set.py` takes seconds.

**Memory.** Use `--workers 4` or fewer for indels: each worker holds the full pooled
candidate set, and the 697 MB sQTL pool OOMs a 125 GB machine at 12 workers. Run
sQTL separately at `--workers 2`.

| Flag | Default | Meaning |
|---|---|---|
| `--pip` | 0.9 | PIP threshold for positives |
| `--exclude_pip` | 0.01 | PIP at which a (variant, gene) pair is excluded from the negative pool and the variant is added to the cross-tissue blocklist |
| `--workers` | 4 | Parallel tissue workers |
| `--qtl` | `eqtl sqtl paqtl` | Modalities to run |
| `--indel_t` | 0 | 0 = SNPs only; >0 = indels only, up to this length |
| `--out_root` | *required* | Root under which `{modality}/` set directories are written |
| `--cache_dir` | `data/gtex11/portal/cache` | Where the two shared preprocessing artifacts live |

### Pipeline overview

`make_vcfs.py` drives the three per-tissue scripts in parallel. For each modality,
two preprocessing artifacts are computed once and reused:

1. **`high_pip_variants.txt`** — every variant whose SuSiE PIP reaches
   `--exclude_pip` for any protein-coding phenotype in any tissue. Used as a
   cross-tissue blocklist, so a weak signal in tissue A cannot be selected as a
   negative in tissue B.
2. **`signif_pool.parquet`** — union of all tissues' `signif_pairs` collapsed to one
   row per (variant, phenotype), with `af` reduced to its median across the tissues
   in which the pair was significant. Pooling gives low-power tissues a denser
   candidate set than their own `signif_pairs` would.

Both depend only on `(modality, --exclude_pip)`, so they live in
`--cache_dir/{qtl}_excl{N}/` and are shared across builds. At ~940 MB they dominate
the pipeline's footprint. They are reused if present; delete to force recomputation.

After preprocessing, per-tissue jobs fan out and their outputs are merged.

### Per-tissue construction

Per tissue, each script:

1. Loads the SuSiE summary parquet, restricts to `biotype == "protein_coding"`,
   parses GTEx-style `variant_id` (`chr_pos_ref_alt_b38`), applies the variant-class
   filter, and computes `maf = min(af, 1 − af)`.
2. Joins per-gene median TPM for the tissue; genes with no TPM here are dropped.
3. Builds modality-specific feature annotation from the GTF.
4. Selects positives at `pip ≥ pos_pip` and annotates each with its feature distance.
5. Loads and filters the candidate negative pool.
6. Greedy-matches each positive to one negative under a feature-aware score.
7. Writes per-tissue VCFs and a `matches.tsv`.

#### Allele matching

SNPs bucket on exact `(ref, alt)`. With only twelve combinations every bucket stays
dense, so exactness is free.

Indels cannot: the allele space is far too large. They bucket instead on
`ilen = len(alt) − len(ref)`, the signed length change, which pins the event type and
size — the property that displaces the model's output bins. Nucleotide identity moves
into the score as `--allele_weight × mismatch`, the Hamming distance between the two
**event alleles**, where the event allele is the longer of REF/ALT: the inserted or
deleted bases plus their anchor base. All members of a bucket share an event-allele
length, so this is well defined.

A CAG insertion can therefore match a CTG insertion for a penalty of 1, but never a
2 bp insertion or any deletion. Exact-allele candidates score 0 and win unless a
mismatched candidate is better by more than `--allele_weight` on the covariate terms.
For SNPs the term is identically zero, so SNP sets are unaffected.

#### eQTL (`src/eqtl_vcfs.py`)

**Region annotation.** Each (variant, gene) pair is labeled against the union of
basic-annotation transcripts of the eQTL gene:

- **UTR3** — variant overlaps any 3′ UTR. 3′ UTR is derived per transcript from CDS
  span + strand: a UTR feature is 3′ when it lies strictly past the CDS in the
  transcribed direction.
- **CDS** — variant overlaps any CDS.
- **TSS** — everything else; carries `tss_dist`, the minimum distance to any TSS of
  the gene.

When a variant overlaps both UTR3 and CDS across different transcripts of the same
gene, the canonical transcript (MANE_Select preferred, GENCODE_Primary fallback)
breaks the tie; if it covers neither, UTR3 wins. A variant mapping to multiple genes
collapses to one row with precedence `UTR3 > CDS > TSS`, ties broken on smallest
`tss_dist`.

**Positives.** SuSiE rows with `pip ≥ pos_pip` and a numeric `afc`; rows without an
AFC estimate are dropped (rare).

**Negatives.** From `signif_pool.parquet`, restricted to protein-coding genes present
in this tissue's SuSiE summary, with focal-tissue positives and blocklisted variants
removed. For **TSS-region negatives only**, the gene assignment is replaced by the
owner of the genome-wide nearest TSS (`pyranges.k_nearest`), and `tss_dist`/`tpm` are
recomputed against it — a candidate's `signif_pairs` gene may have a far-away TSS
while the variant actually sits at a different gene's promoter, and the reassignment
makes the matching distance meaningful.

**Score.** Candidates bucket by `(region, ref, alt)` for SNPs, `(region, ilen)` for
indels, and are ranked by

```
score = |Δ log1p(TPM)| + λ_d · |Δ log1p(tss_dist)| + λ_m · |Δ log(MAF)| + λ_a · mismatch
```

with all weights defaulting to 1. The distance term applies only to `region == TSS`;
in UTR3 and CDS the bucketing already pins the variant to the same compartment.
Lower is better. Greedy assignment: each positive claims the lowest-scoring available
candidate, each negative is used at most once, ties broken stably.

#### sQTL (`src/sqtl_vcfs.py`)

**Splice sites** are the internal exon boundaries of the no-readthrough
protein-coding transcripts: every non-terminal exon end (donor) and every non-first
exon start (acceptor), deduplicated by `(chrom, pos, gene_id)`.

Positives take the gene-specific nearest splice site; positives whose gene has no
splice site are dropped, and junction-level duplicates collapse to maximum PIP.
Negatives take their **genome-wide nearest** splice site, with the gene reassigned to
that site's owner for the TPM lookup, deduplicated to one row per variant by minimum
distance. Buckets are `(ref, alt)` / `(ilen,)`; the score has the same form as eQTL
with `splice_dist` in place of `tss_dist`.

#### apaQTL (`src/paqtl_vcfs.py`)

Parallel to sQTL with PAS sites — the annotated 3′ end of each transcript in the
no-readthrough GTF, gene-tagged — replacing splice sites. Bucketing and score
identical in form.

### Reproducibility

- `np.random.seed(0)` is set inside each per-tissue script.
- `argsort(..., kind="stable")` and `groupby(..., sort=False)` keep tie-breaking and
  bucket ordering deterministic given fixed input parquets.
- Per-tissue jobs are independent processes, so worker count does not affect outputs.
- The GTF is rebuilt from GENCODE by `portal.mk` rather than shipped, so the
  annotation provenance is in this repo.
- Positive VCFs are sorted by PIP descending, ties broken by `variant_id`, so row
  order is canonical rather than incidental. Consumers should still key on the ID
  column: `{tissue}_pos.vcf` and a score file computed from it correspond by row,
  and pairing a rebuilt VCF with a stale score file is not detected.

## QC

[src/qc.ipynb](src/qc.ipynb) covers all three modalities side by side.
**Set-level** (`_pos.vcf`): counts, MAF/TPM/distance/effect/PIP distributions, MAD
outlier heatmap, KS tests, anomalous-tissue ranking. **Pair-level**
(`matches.tsv`): match rate, pos-vs-neg marginals and scatter, same-gene rate, score
by region and tissue — marginals can agree while the pairing is scrambled, which only
the scatter reveals. [src/qc_indels.ipynb](src/qc_indels.ipynb) runs the same checks
for the indel sets, each statistic beside its SNP counterpart.

Both are generated by `src/mk_qc.py`; edit that, not the `.ipynb`.

### Known limitation: low-powered tissues

GTEx tissue sample size, the number of high-PIP fine-mapped variants, and the
biological plausibility of those variants move together. Tissues with small cohorts —
primarily kidney, bladder, vagina, and uterus — produce positive sets that look
anomalous relative to well-powered tissues.

The clearest case is **bladder** (121 eQTL positives vs ~2,700 for Nerve_Tibial). Its
positives have a median TSS distance of ~153 kb versus ~9 kb for Nerve_Tibial, and a
meaningful fraction sit more than 500 kb from any TSS of their target gene. Those
very distant bladder positives are enriched for low minor allele frequencies
(Mann-Whitney and KS both p < 0.05 against nearer variants in the same tissue). The
pattern is consistent with a statistical artifact: with few samples, a handful of
individuals sharing an allele *and* high expression can drive SuSiE to high PIP
without the variant having any mechanistic role near the gene.

For cross-tissue evaluation the effect is diluted, since the affected tissues
contribute a small fraction of the pool. For tissue-specific work, prefer the
[gold sets](#gold-sets), which address both halves of this: the distance cap removes
the far-TSS artifacts, and the tissue floor drops bladder and kidney cortex in all
three modalities.

## Releases and versioning

Three things are named, one each:

| layer | name | what it is |
|---|---|---|
| repo / method | `qtl-bench` | the construction pipeline, gold filter, and QC |
| study / cohort | `gtex11` | the input data source |
| release | `gtex11_v1` | one build of the study by the method |

A release is a git tag, a `{tag}.tar.gz` asset, and a
[CHANGELOG.md](CHANGELOG.md) section. There are no version directories: the live sets
sit at `data/{study}/`, and `VERSION` names the release they came from. To compare
against an older release, unpack it somewhere else —
`./fetch.sh gtex11_v1 data/gtex11/archive`.

Set directory names carry identity only (study, variant class, modality, tier);
`build.json` and `gold.json` carry the parameters. A different PIP threshold is a new
release, not a new directory name.

## Repository layout

```
fetch.sh          download a release
portal.mk         download every public input
src/
  make_vcfs.py    orchestrator: preprocessing, parallel fan-out, merge
  sc/             single-cell studies: fetchers, input preparation, builders (qb_sets.py
                  shared matching)
  {eqtl,sqtl,paqtl}_vcfs.py   per-tissue construction
  gold_set.py     derive a _gold subset of a set directory
  qc_lib.py       VCF / match-table loaders shared by the notebooks
  mk_qc.py        generates qc.ipynb
  qc.ipynb  qc_indels.ipynb
docs/
  methods.md          paper-style Methods section
  gold_criteria.md    evidence behind the gold thresholds
data/             gitignored; releases unpack here
```

## Citation

If you use these sets, please cite this repository and the release the study is built
on: the GTEx Consortium v11 release (`gtex11`); the ADSP Functional Genomics xQTL Atlas,
NIAGADS NG00184 (`adsp_brain`); or TenK10K Phase 1 (`tenk10k_pbmc`): Cuomo et al. for the
eQTLs, and Xue et al., medRxiv 2025.08.27.25334533, for the caQTLs. The underlying summary
statistics and fine-mapping are theirs; this repository contributes the positive
selection, negative matching, and gold filtering.

## License

Code is Apache-2.0 ([LICENSE](LICENSE)). The released benchmark data is CC-BY-4.0,
derived from open-access summary statistics: GTEx v11, NIAGADS NG00184 (open access;
see its data use terms), and TenK10K Phase 1 (CC-BY-4.0). The TenK10K caQTL release
accompanies a preprint that had not completed peer review when these sets were built.
