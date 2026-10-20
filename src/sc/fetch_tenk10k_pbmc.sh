#!/bin/bash
# Fetch the TenK10K inputs: eQTL nominal cis statistics (the signif_pairs analogue qtl-bench
# draws negatives from) and caQTL fine-mapping. Resumable: re-run to continue.
#
#   eQTL:   Zenodo 10.5281/zenodo.18221260, all cis tests at 100 kb (16.6 GB) plus the
#           gene-level results (per-gene thresholds, 10 MB)
#   caQTL:  HF anglixue/TenK10K_multiome, SuSiE fine-mapping (positives only; no negatives)
#
#   fetch_tenk10k_pbmc.sh [eqtl|caqtl|all] [portal dir]   (default data/tenk10k_pbmc/portal;
#   the two run fine in parallel)
set -euo pipefail

what=${1:-all}
P=${2:-$(dirname $0)/../../data/tenk10k_pbmc/portal}

if [ $what = eqtl ] || [ $what = all ]; then
Z=https://zenodo.org/records/18221260/files
mkdir -p $P/eqtl/raw
cd $P/eqtl
# Zenodo answers a resume of a complete file with 416, so skip those by size.
while read f size; do
    [ "$(stat -c %s $f 2>/dev/null)" = $size ] && continue
    curl -fL -C - --retry 5 -o $f "$Z/$f?download=1"
done <<'EOF'
common_variant_gene_level_results_with_annotated_fails.zip 10137272
common_all_cis_pvalues_100kb.zip 16559137990
EOF
curl -s https://zenodo.org/api/records/18221260 > raw/zenodo_record.json
fi

if [ $what = caqtl ] || [ $what = all ]; then
HF=https://huggingface.co/datasets/anglixue/TenK10K_multiome/resolve/main
mkdir -p $P/caqtl/raw
cd $P/caqtl
[ -s all_celltypes_caQTL_susie_annotated.zip ] || \
    curl -fL --retry 5 -o all_celltypes_caQTL_susie_annotated.zip "$HF/Fine_mapping/all_celltypes_caQTL_susie_annotated.zip"
echo "b11121fc6300df42afd78fa1c23efdd4  all_celltypes_caQTL_susie_annotated.zip" | md5sum -c -
unzip -oq all_celltypes_caQTL_susie_annotated.zip -x '__MACOSX/*' -d raw/
fi
echo done
