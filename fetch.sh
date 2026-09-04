#!/bin/sh
# Download a benchmark release and unpack it.
#
#   ./fetch.sh gtex11_v1                        -> data/gtex11/
#   ./fetch.sh gtex11_v1 data/gtex11/archive    -> somewhere else, to compare
set -e

tag=$1
if [ -z "$tag" ]; then
    echo "usage: $0 <release-tag> [dest]" >&2
    exit 1
fi
dest=${2:-data/${tag%_v*}}

mkdir -p "$dest"
curl -fL "https://github.com/calico/qtl-bench/releases/download/$tag/$tag.tar.gz" \
    | tar xz -C "$dest"
echo "unpacked $tag into $dest"
