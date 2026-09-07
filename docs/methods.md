# Methods: Construction of fine-mapped GTEx variant benchmark sets

## Overview

To evaluate sequence-to-function models on cis-regulatory variant effect prediction, we constructed positive and matched-negative variant sets from the GTEx v11 release for three molecular phenotypes — gene expression (eQTL), pre-mRNA splicing (sQTL), and alternative polyadenylation (apaQTL) — across all 50 GTEx v11 tissues. Positives are statistically fine-mapped causal variants (SuSiE posterior inclusion probability ≥ 0.9); negatives are drawn from the same modality's pool of nominally-significant cis variants and are matched one-to-one to each positive on local genomic context, expression of the linked gene, distance to the nearest cognate regulatory feature, minor allele frequency, and exact reference and alternate alleles. Construction runs from a single driver script (`src/make_vcfs.py`) and is deterministic: given the same GTEx input parquets and annotation, a rebuild reproduces the same variants, pairings, and scores.

## Data sources

We used GTEx v11 cis-QTL summary statistics and SuSiE fine-mapping released by the GTEx Consortium (https://www.gtexportal.org). For each of the three phenotype modalities (eQTL, sQTL — LeafCutter intron-cluster usage; apaQTL — 3′ UTR usage ratios), we downloaded both the `signif_pairs` parquet files (nominally-significant variant–phenotype pairs from tensorQTL) and the `SuSiE_summary` parquet files (per-(variant, phenotype) PIP, credible-set membership, and where applicable allelic fold change) for all 50 tissues. We also used the GTEx v11 RNA-SeQC v2.4.3 gene-level median-TPM matrix.

Gene annotation used GENCODE v48 basic protein-coding transcripts (no-readthrough subset) on GRCh38. Excluding readthrough transcripts avoids contaminating UTR3 annotations with spuriously extended 3′ UTRs, and keeps splice-site and polyadenylation-site sets clean of readthrough-specific junctions and 3′ ends. Canonical-transcript tags were taken directly from the GTF: `MANE_Select` was preferred and `GENCODE_Primary` was used as a fallback for genes that lack a MANE assignment.

## Variant filters

Across all modalities, we restricted the SuSiE summary to phenotypes assigned `biotype == "protein_coding"`, dropped variants with missing median TPM in the focal tissue, and retained SNPs only. Minor allele frequency was computed as `MAF = min(AF, 1 − AF)` from the GTEx v11 alternate-allele frequency.

Insertions and deletions were built into a separate, parallel set of benchmarks using an otherwise identical construction. Indels were capped at 4 bp, since longer indels displace the model's output bins enough to complicate scoring. All GTEx v11 indels are left-anchored single-event alleles (no complex substitutions or MNPs), so no allele normalization was required. Indels are bucketed by signed length change and charged a per-base allele mismatch penalty rather than bucketed on exact alleles (see *Score-based one-to-one matching*). Region overlap for eQTL was computed over the full reference-allele span rather than a single base, so multi-base deletions straddling a compartment boundary are placed correctly.

### Gold subsets

For faster evaluation on unambiguous positives, we additionally derived a *gold* subset of each set by filtering positives on three criteria and dropping each rejected positive's matched negative with it. First, distance: positives more than 500 kb (eQTL) or 100 kb (sQTL, apaQTL) from the cognate feature were removed, as they fall outside the receptive field of the models being evaluated; eQTL positives inside the gene body carry no TSS distance and were retained unconditionally. Second, significance: positives with no nominal p-value in `signif_pairs` were removed, since they are SuSiE secondary signals detectable only after conditioning on the primary, and are therefore drawn from a different population than the negatives they are matched to; the remainder were required to reach −log10 p ≥ 8. Third, confidence: PIP = 1.0, which restricts to credible sets that SuSiE resolved to a single variant. Finally, tissues retaining fewer than one-third of the modality's median surviving positive count were dropped whole, which removes the low-powered cohorts (bladder, kidney cortex, amygdala) discussed above.

## Positive variant definition

For each tissue and modality, a positive variant–phenotype pair was defined as a row in the SuSiE summary with PIP ≥ 0.9. For eQTL we additionally required a numeric SuSiE allelic fold change; rows without an AFC estimate were dropped. Per-(tissue, gene) nominal p-values and per-allele slopes (sQTL/apaQTL) were looked up from the corresponding `signif_pairs` parquet for inclusion in the output VCF.

When a variant entered the positive set against multiple phenotypes — for example, overlapping eQTL credible sets for adjacent protein-coding genes in the same cis window — we collapsed to one row per variant. For eQTL this collapse uses region precedence (3′ UTR > CDS > TSS, see below); for sQTL and apaQTL it keeps the maximum-PIP row.

## Local context annotation

**eQTL.** Each (variant, gene) positive was labeled as `UTR3`, `CDS`, or `TSS`-proximal relative to its eQTL gene:
- *UTR3* if the variant lies inside any 3′ UTR interval of the gene. 3′ UTRs were derived per transcript from CDS span and strand.
- *CDS* if the variant lies inside any CDS interval of the gene.
- *TSS* otherwise; carries `tss_dist`, the minimum absolute distance to any TSS of the gene.

When a variant overlapped both a 3′ UTR and a CDS interval across different transcripts, the canonical transcript (MANE_Select preferred, GENCODE_Primary fallback) determined the label; if the canonical transcript covered neither, the variant was assigned UTR3 by default. Cross-gene multi-mappings collapsed to a single row with precedence UTR3 > CDS > TSS, ties broken by smallest TSS distance.

**sQTL.** Splice sites were defined as the internal exon boundaries of the no-readthrough basic protein-coding transcripts: every non-terminal exon end (donor) and every non-first exon start (acceptor), per transcript, deduplicated by `(chrom, position, gene)`. Each positive was annotated with the minimum distance to a splice site of its gene; positives whose gene had no splice site were dropped.

**apaQTL.** Polyadenylation sites were proxied by the annotated 3′ end of every transcript in the no-readthrough GTF, gene-tagged. Each positive was annotated with the minimum distance to a PAS of its gene.

## Negative candidate pool

Per modality, we built a single cross-tissue candidate pool by unioning all 50 tissues' `signif_pairs`, deduplicating to one row per `(variant_id, phenotype_id)`, and taking the median allele frequency across tissues in which the pair was observed. Pooling provides a denser candidate set for low-power tissues that would otherwise have sparse `(region, ref, alt)` buckets.

For each focal tissue the pool was restricted to protein-coding genes with nonzero median TPM and at least one SuSiE row in that tissue. We excluded (a) focal-tissue positives and (b) variants present in the cross-tissue blocklist — every variant with SuSiE PIP ≥ 0.01 for any protein-coding phenotype in any tissue, computed in a one-pass scan over all SuSiE summaries before per-tissue jobs begin. The cross-tissue blocklist prevents weakly fine-mapped variants in tissue A from serving as negatives in tissue B.

## Negative gene and distance reassignment

A negative candidate's `signif_pairs` gene may have a regulatory feature far from the variant. To ensure the matching distance reflects the variant's nearest plausibly-regulated gene, we reassigned each negative's gene to the owner of the genome-wide nearest cognate feature (via `pyranges.k_nearest`): nearest TSS for eQTL TSS-region candidates; nearest splice site for sQTL; nearest PAS for apaQTL. TPM was reloaded from the reassigned gene. UTR3 and CDS eQTL negatives keep their original gene assignment because compartment overlap, not distance, drives their matching.

## Score-based one-to-one matching

Within each tissue and modality, we paired each positive to one negative by first bucketing on alleles and (for eQTL) region. For SNPs the allele key is exact — `(region, ref, alt)` for eQTL, `(ref, alt)` for sQTL and apaQTL — which is free to enforce because there are only twelve reference/alternate combinations. For indels the allele space is far too large for exact bucketing to fill every bucket, so the key is the signed length change `Δ = len(alt) − len(ref)`, which fixes the event type and size, and nucleotide identity is charged to the score instead.

Inside each bucket, candidates were ranked by

```
s = |Δ log1p(TPM)| + λ_d · |Δ log1p(d)| + λ_m · |Δ log(MAF)| + λ_a · h
```

where `d` is the modality-specific feature distance, `λ_d = 1` is the distance weight, `λ_m = 1` is the MAF weight, and `h` is the Hamming distance between the positive's and candidate's event alleles under `λ_a = 1`. The event allele is the longer of REF and ALT — the inserted or deleted bases together with their anchor base — and is of equal length for all members of a bucket. A three-base insertion of CAG can therefore be matched to an insertion of CTG at a cost of one, but never to a two-base insertion or to a deletion. Exact-allele candidates score `h = 0` and are preferred unless a mismatched candidate is better by more than `λ_a` in the covariate terms; for SNPs, whose buckets are already allele-exact, `h` is identically zero. The distance term was suppressed for UTR3 and CDS eQTL buckets, where compartment overlap already controls for genomic context. Lower scores are better. Greedy assignment: each positive claims the lowest-scoring available candidate in its bucket, with each negative used at most once.

## Output

Per modality and tissue: `{tissue}_pos.vcf`, `{tissue}_neg.vcf`, and `{tissue}_matches.tsv`. Across all tissues: a deduplicated sorted union VCF (`merge.vcf`) and a concatenated matches table; the separate `pos_merge.vcf` / `neg_merge.vcf` unions are retained for backward compatibility but are deprecated. The PIP ≥ 0.9 SNP release contains 58,122 / 58,517 / 14,823 tissue-level pairs and 23,836 / 20,732 / 4,991 unique positive variants for eQTL / sQTL / apaQTL; the gold subsets contain 16,415 / 19,371 / 4,548 pairs across 42 / 47 / 47 tissues. VCF INFO fields encode gene, tissue, PIP, region (eQTL) or feature distance (sQTL/apaQTL), TPM, MAF, allelic fold change (eQTL positives) or per-allele slope (sQTL/apaQTL positives), and the −log10 nominal p-value.
