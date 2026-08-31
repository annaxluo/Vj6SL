#!/usr/bin/env bash
set -euo pipefail

## Install gffread 0.12.7 
##
## Usage:
##   export SOFTWARE_ROOT=/path/to/shared/software
##   bash modules/gffread/install_gffread_0.12.7.sh

VERSION="0.12.7"
TOOL="gffread"

SOFTWARE_ROOT="${SOFTWARE_ROOT}"
INSTALL_DIR="${SOFTWARE_ROOT}/${TOOL}/${VERSION}"
BUILD_DIR="${TMPDIR:-/tmp}/${TOOL}-${VERSION}-build-${USER:-user}-$$"

mkdir -p "${INSTALL_DIR}/bin"
mkdir -p "${BUILD_DIR}"

echo "Installing ${TOOL} ${VERSION}"
echo "SOFTWARE_ROOT: ${SOFTWARE_ROOT}"
echo "INSTALL_DIR:   ${INSTALL_DIR}"
echo "BUILD_DIR:     ${BUILD_DIR}"

cd "${BUILD_DIR}"

wget "https://github.com/gpertea/gffread/releases/download/v${VERSION}/gffread-${VERSION}.Linux_x86_64.tar.gz" \
    -O "gffread-${VERSION}.Linux_x86_64.tar.gz"

tar -xzf "gffread-${VERSION}.Linux_x86_64.tar.gz"

if [[ -f "gffread-${VERSION}.Linux_x86_64/gffread" ]]; then
    cp "gffread-${VERSION}.Linux_x86_64/gffread" "${INSTALL_DIR}/bin/"
elif [[ -f "gffread" ]]; then
    cp "gffread" "${INSTALL_DIR}/bin/"
else
    echo "Error: could not find gffread binary after extraction."
    exit 1
fi

chmod -R a+rX "${INSTALL_DIR}"
chmod a+rx "${INSTALL_DIR}/bin/gffread"

echo "gffread installed to: "
echo "${INSTALL_DIR}"
