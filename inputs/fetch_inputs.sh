#!/usr/bin/env bash
# Obtain the exact one- and two-loop Type-II finite remainders, published as ancillary
# files of arXiv:1909.08993, and check them against inputs/SHA256SUMS.
#
#   bash inputs/fetch_inputs.sh                     # download from arXiv
#   bash inputs/fetch_inputs.sh path/to/source.tgz  # use a local copy of the arXiv source
set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
target="$here/arXiv-1909.08993"
mkdir -p "$target"

if [[ $# -ge 1 ]]; then
  tarball=$1
else
  tarball="$target/source.tar.gz"
  curl --fail --location --output "$tarball" https://arxiv.org/src/1909.08993v4
fi

tar -xzf "$tarball" -C "$target" anc
(cd "$target/anc" && sha256sum --check "$here/SHA256SUMS")
echo "Inputs ready in $target/anc"
