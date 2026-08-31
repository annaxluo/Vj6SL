#!/usr/bin/env bash
set -euo pipefail

# Install STAR 2.7.8a.
#
# Usage:
#   export SOFTWARE_ROOT=/path/to/shared/software
#   bash modules/star/install_star_2.7.8a.sh

VERSION="2.7.8a"
INSTALL_DIR="${SOFTWARE_ROOT}/STAR/${VERSION}"
BUILD_DIR="${INSTALL_DIR}/build"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

wget https://github.com/alexdobin/STAR/archive/${VERSION}.tar.gz
tar -xzf ${VERSION}.tar.gz

cd "STAR-${VERSION}/source"
make

mkdir -p "${INSTALL_DIR}/bin"
cp STAR "${INSTALL_DIR}/bin/STAR"

chmod -R a+rX "${INSTALL_DIR}"

echo "STAR ${VERSION} installed in:"
echo "${INSTALL_DIR}"
echo "STAR binary:"
echo "${INSTALL_DIR}/bin/STAR"
