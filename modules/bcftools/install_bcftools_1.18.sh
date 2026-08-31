#!/usr/bin/env bash
set -euo pipefail

## Install bcftools 1.18 
##
## Usage:
##   export SOFTWARE_ROOT=/path/to/shared/software
##   bash modules/bcftools/install_bcftools_1.18.sh

VERSION="1.18"
TOOL="bcftools"

SOFTWARE_ROOT="${SOFTWARE_ROOT}"
INSTALL_DIR="${SOFTWARE_ROOT}/${TOOL}/${VERSION}"
BUILD_DIR="${TMPDIR:-/tmp}/${TOOL}-${VERSION}-build-${USER:-user}-$$"

mkdir -p "${INSTALL_DIR}"
mkdir -p "${BUILD_DIR}"

cd "${BUILD_DIR}"

wget "https://github.com/samtools/bcftools/releases/download/${VERSION}/bcftools-${VERSION}.tar.bz2" \
    -O "bcftools-${VERSION}.tar.bz2"

tar -xjf "bcftools-${VERSION}.tar.bz2"
cd "bcftools-${VERSION}"

./configure --prefix="${INSTALL_DIR}"
make
make install

chmod -R a+rX "${INSTALL_DIR}"

echo "bcftools installed to: "
echo "${INSTALL_DIR}"
