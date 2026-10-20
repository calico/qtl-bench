#!/bin/bash
# Fetch ADSP FunGen xQTL Atlas (NIAGADS NG00184.v1, open access) snuc-eQTL archives needed
# to replicate qtl-bench negatives: significance-filtered pairs (the signif_pairs analogue)
# and all-variant fine-mapping (the PIP >= 0.01 blocklist needs PIPs outside cs95).
# File IDs from https://st1.niagads.org/portal/v1/fileset_public/NG00184. Resumable.
set -euo pipefail

# fetch_adsp_brain.sh [portal dir]   (default data/adsp_brain/portal)
out=${1:-$(dirname $0)/../../data/adsp_brain/portal}/raw/NG00184
mkdir -p $out
cd $out
API=https://st1.niagads.org/portal/v1/download-public/NG00184/fileset
curl -sL -o NG00184.v1.file_manifest.csv https://st1.niagads.org/portal/download-public/NG00184.v1/fm

# file_ID name md5
while read id name md5; do
    curl -fL -C - --retry 5 -o $name $API/$id
    [ "$md5" = - ] || echo "$md5  $name" | md5sum -c -
done <<'EOF'
57106 ADSP_FunGen_xQTL.v1.snuc-eQTL.hmt_significant.manifest.txt 2e8e80a8a0c8cac321688a26c7789b15
57104 ADSP_FunGen_xQTL.v1.snuc-eQTL.bh_significant.manifest.txt 176e232fe99a2fb3f0cfcd72dc4ecc05
57108 ADSP_FunGen_xQTL.v1.snuc-eQTL.single_context_finemapping_all.manifest.txt 7f6be77e94eca9ed825ab5f2498d5af3
57107 ADSP_FunGen_xQTL.v1.snuc-eQTL.hmt_significant.tar fa0008baccf2067217fb84d8a378c27d
57109 ADSP_FunGen_xQTL.v1.snuc-eQTL.single_context_finemapping_all.tar f0d3457b1ed556f85cfbc2651a20190f
57105 ADSP_FunGen_xQTL.v1.snuc-eQTL.bh_significant.tar fed69bc64e7cbd888bedaf06261b4dba
EOF
echo done
