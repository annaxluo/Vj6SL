#!/usr/bin/env bash
set -euo pipefail

# htslib 1.18 installation script.
#
# Usage:
#   export SOFTWARE_ROOT=/path/to/shared/software
#   bash modules/htslib/install_htslib_1.18.sh

VERSION="1.18"
TOOL="htslib"

SOFTWARE_ROOT="${SOFTWARE_ROOT}"
INSTALL_DIR="${SOFTWARE_ROOT}/${TOOL}/${VERSION}"
BUILD_DIR="${TMPDIR:-/tmp}/${TOOL}-${VERSION}-build-${USER:-user}-$$"

mkdir -p "${INSTALL_DIR}"
mkdir -p "${BUILD_DIR}"

cd "${BUILD_DIR}"

wget "https://github.com/samtools/htslib/releases/download/${VERSION}/htslib-${VERSION}.tar.bz2" \
    -O "htslib-${VERSION}.tar.bz2"

tar -xjf "htslib-${VERSION}.tar.bz2"
cd "htslib-${VERSION}"

./configure --prefix="${INSTALL_DIR}"
make
make install

chmod -R a+rX "${INSTALL_DIR}"

echo "htslib installed to: "
echo "${INSTALL_DIR}"

