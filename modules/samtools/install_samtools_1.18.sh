#!/usr/bin/env bash
set -euo pipefail

# Install Samtools 1.18.
#
# Usage:
#   export SOFTWARE_ROOT=/path/to/shared/software
#   bash modules/samtools/install_samtools_1.18.sh

VERSION="1.18"
INSTALL_DIR="${SOFTWARE_ROOT}/samtools/${VERSION}"
BUILD_DIR="${INSTALL_DIR}/build"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

wget https://github.com/samtools/samtools/releases/download/${VERSION}/samtools-${VERSION}.tar.bz2 -O samtools.tar.bz2
tar -xjf samtools.tar.bz2
rm samtools.tar.bz2

cd "samtools-${VERSION}"

./configure --prefix="${INSTALL_DIR}"
make
make install

cd "${INSTALL_DIR}"

rm -rf "${BUILD_DIR}"

chmod -R a+rX "${INSTALL_DIR}"

echo "Samtools ${VERSION} installed in:"
echo "${INSTALL_DIR}"