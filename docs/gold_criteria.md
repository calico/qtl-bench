# Gold-set criteria

`src/gold_set.py` filters a benchmark set down to positives whose causal identity and
measured effect are both unambiguous, so that model evaluation is faster and its
failures are attributable to the model rather than to the label. This document
records why each threshold is where it is.

| criterion | eQTL | sQTL | apaQTL |
|---|---|---|---|
| distance cap | `TSSD` < 500 kb | `SD` < 100 kb | `PD` < 100 kb |
| significance | `NLP` present and ≥ 8 | same | same |
| confidence | `PIP` = 1.0 | same | same |
| tissue floor | ≥ 0.33 × the modality's median surviving count | same | same |

Realized:

| | full set | gold | reduction | median positives/tissue |
|---|---|---|---|---|
| eQTL | 58,122 · 50 tissues | 16,415 · 42 | 3.5× | 1,140 → 350 |
| sQTL | 58,517 · 50 | 19,371 · 47 | 3.0× | 1,084 → 350 |
| apaQTL | 14,823 · 50 | 4,548 · 47 | 3.3× | 262 → 82 |

**Each criterion is structural first.** Choosing a benchmark by one model's accuracy
is circular, so every cut below has a justification that does not appeal to any
model: >500 kb is outside Borzoi's 524 kb receptive field; a positive with no
nominal p-value is not drawn from the pool its negative is sampled from; `PIP < 0.95`
means SuSiE's credible set holds more than one variant. The AUPRC tables are
corroboration. The one criterion resting mainly on measurement is `PIP = 1.0` over
`PIP ≥ 0.95`; that is stated plainly below.

Measurements are an internal Borzoi-class ensemble scored with `covgene-logFC` and
`covgene-lnDinf` over 48 eval tissues, AUPRC of positives against their matched
negatives. Counts in the AUPRC tables are the scored eval pool, slightly smaller
than the full set.

## Distance

Beyond the caps the model is at chance:

| distance bin | eQTL | sQTL | apaQTL |
|---|---|---|---|
| < 1 kb | 0.801 | 0.747 | 0.758 |
| 1–10 kb | 0.763 | 0.766 | 0.726 |
| 10–100 kb | 0.742 | 0.762 | 0.738 |
| 100–500 kb | 0.598 | **0.505** | **0.482** |
| > 500 kb | **0.498** | — | — |

The eQTL cap stays at 500 kb rather than tightening to 100 kb: it removes only the
unscoreable bin and keeps the long-range enhancer biology a benchmark should still
ask about, even though 100–500 kb is already degraded at 0.598.

eQTL exonic positives (`REGION` in `CDS`/`UTR3`) carry no `TSSD` — they are inside
the gene by construction — and pass unconditionally.

## Missing `NLP`

`NLP` is a positives-only INFO field; every negative row carries `NLP=.` by design,
so the asymmetry is not visible in the VCF. It sits one level down. Negatives are
*sampled from* `signif_pool.parquet`, the union of `signif_pairs`, so every negative
has a nominal association by construction. A positive with no `NLP` is absent from
`signif_pairs` entirely, so the two members of the pair are drawn from different
populations, which breaks the premise of the matched design.

This is not leakage in the existing sets, which was checked directly: all 1,308
missing-`NLP` eQTL positives appear in `high_pip_variants.txt` (`collect_high_pip`
scans SuSiE summaries at `pip >= 0.01`, catching them whether or not they have a
p-value), and 0 of 46,754 unique negatives appear in the positive set or the
blocklist. The question is only which positives belong in a gold set.

These positives are SuSiE *secondary* signals (`cs_id ≥ 2`), detectable only after
conditioning out the primary, which a marginal tensorQTL scan never sees. They are
genuine eQTLs and are retained in the full sets, but a model scored on marginal
effect cannot be expected to resolve a conditional-only signal.

The measurement agrees: 1,779 / 1,472 / 492 positives (3.3% / 2.7% / 3.5%) at AUPRC
0.636 / 0.602 / 0.627 against 0.702 / 0.683 / 0.674 for the rest.

## `NLP ≥ 8`

`signif_pairs` is already p-thresholded, so there is essentially nothing below
`NLP` 4 and the cut removes the bottom two bins. Within the distance cap, the step
from the `6–8` bin to `8–12` is clear in two modalities (eQTL 0.660 → 0.696, apaQTL
0.649 → 0.690) and gradual in sQTL (0.646 → 0.665). Going on to `≥15` costs ~6,000
eQTL positives for +0.007. 8 is the knee, and it is the `NLP_WEAK` constant
[qc.ipynb](../src/qc.ipynb) already flags on.

`NLP` is not redundant with `PIP`: on positives they correlate only weakly
(Spearman 0.34 / 0.42 / 0.38), and within `PIP ≥ 0.999` the `NLP` ladder still
climbs — present 0.746, `≥8` 0.750, `≥15` 0.757, `≥25` 0.765.

## `PIP = 1.0`

**The structural floor is 0.95, not 0.9.** Across 8 tissues' SuSiE summaries
(1.17 M protein-coding rows), every row with `PIP ≥ 0.95` has `cs_size == 1`, while
the 0.90–0.95 band is 0% singleton with a mean credible set of 2.5 variants. Those
are variants SuSiE could not distinguish from an LD partner — disqualifying on
principle for a set whose premise is "this variant is the causal one," whatever the
AUPRC says.

**The empirical knee is at exactly 1.0**, not 0.999 (within the distance cap, `NLP`
present):

| PIP | eQTL | sQTL | apaQTL |
|---|---|---|---|
| .90–.95 | 0.670 | 0.676 | 0.642 |
| .95–.99 | 0.672 | 0.679 | 0.657 |
| .99–.999 | 0.700 | 0.691 | 0.674 |
| .999–.9999 | 0.716 | 0.704 | 0.706 |
| .9999–<1 | 0.720 | 0.721 | 0.711 |
| **= 1.0** | **0.757** | **0.751** | **0.751** |

`PIP ≥ 0.999` is 70% `PIP = 1.0` diluted by a ~0.72 tail, which is why it lands
between the two. `= 1.0` is both the sharper cut and the more legible criterion:
SuSiE placed the entire posterior on one variant. (VCF `PIP` is written to four
decimals, so the criterion is `PIP ≥ 0.99995` in the underlying fine-mapping.)

**It is real information, not an effect-size proxy.** Spearman(`PIP`, `|AFC|`) =
0.197, and within every `|AFC|` tercile the `PIP` ladder gains a similar ~3.5
points:

| `|AFC|` tercile | PIP .90–.99 | PIP = 1.0 |
|---|---|---|
| low | 0.679 | 0.715 |
| mid | 0.732 | 0.763 |
| high | 0.753 | 0.775 |

**No `|AFC|` floor**, even though `|AFC|` predicts well on its own (0.672 below 0.25,
rising to 0.751 above 2): filtering on measured effect magnitude makes the benchmark
easier by construction rather than cleaner. `PIP` and `NLP` are confidence measures;
`|AFC|` is the quantity being predicted.

## Stacked effect

| | eQTL | sQTL | apaQTL |
|---|---|---|---|
| all | 54,406 · 0.700 | 55,104 · 0.681 | 13,885 · 0.672 |
| + distance | 52,313 · 0.708 | 47,174 · 0.711 | 11,954 · 0.700 |
| + NLP present, ≥ 8 | 41,175 · 0.724 | 41,657 · 0.720 | 10,099 · 0.709 |
| + PIP ≥ 0.99, + tissues | 27,483 · 0.741 | 29,533 · 0.732 | 7,069 · 0.729 |
| + PIP ≥ 0.999, + tissues | 20,860 · 0.750 | 23,390 · 0.742 | 5,542 · 0.741 |
| **+ PIP = 1.0, + tissues** | **15,377 · 0.759** | **18,248 · 0.751** | **4,222 · 0.752** |

A 3–3.5× reduction at +0.06 AUPRC. The one cost is apaQTL thinning to a median of 82
positives per tissue, workable pooled but noisy per tissue; `--pip 0.999` is the
fallback if per-tissue apaQTL numbers matter.

## Tissue floor

An absolute count floor cannot work across modalities — after the variant filters
the per-tissue median is 328 (eQTL), 340 (sQTL) but only 71 (apaQTL) — so the floor
is a fraction of that modality's own median, `--tissue_frac` (default 0.33).
Surviving tissues and the share of positives kept:

| frac | eQTL | sQTL | apaQTL |
|---|---|---|---|
| 0.20 | 46 tissues · 99.3% pos | 49 · 99.7% | 49 · 99.8% |
| **0.33** | **42 · 97.3%** | **47 · 98.6%** | **47 · 99.0%** |
| 0.50 | 35 · 91.8% | 41 · 94.0% | 42 · 95.8% |

**This is hygiene, not a quality or speed lever.** Pooled AUPRC moves by ≤0.001 at
any of these fractions (sQTL 0.7326 → 0.7321, slightly *down*), because the dropped
tissues are too small a share to matter. What it buys is that no surviving tissue is
too thin for a stable *per-tissue* metric. 0.50 overshoots, taking 8.2% of eQTL
positives including Liver.

At 0.33 the realized drops are exactly the tissues the README QC note flags:

| | dropped |
|---|---|
| eQTL (8) | Bladder, Brain_Amygdala, Brain_Substantia_nigra, Kidney_Cortex, Minor_Salivary_Gland, Ovary, Uterus, Vagina |
| sQTL (3) | Bladder, Brain_Amygdala, Kidney_Cortex |
| apaQTL (3) | Bladder, Brain_Amygdala, Kidney_Cortex |
