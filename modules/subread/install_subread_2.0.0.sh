#!/usr/bin/env bash
set -euo pipefail

## Subread 2.0.0 installation.
##
## Usage:
##   export SOFTWARE_ROOT=/path/to/shared/software
##   bash modules/subread/install_subread_2.0.0.sh

VERSION="2.0.0"
TOOL="subread"

SOFTWARE_ROOT="${SOFTWARE_ROOT}"
INSTALL_DIR="${SOFTWARE_ROOT}/${TOOL}/${VERSION}"
BUILD_DIR="${TMPDIR:-/tmp}/${TOOL}-${VERSION}-build-${USER:-user}-$$"

mkdir -p "${INSTALL_DIR}"
mkdir -p "${BUILD_DIR}"

echo "Installing Subread ${VERSION}"
echo "SOFTWARE_ROOT: ${SOFTWARE_ROOT}"
echo "INSTALL_DIR:   ${INSTALL_DIR}"
echo "BUILD_DIR:     ${BUILD_DIR}"

cd "${BUILD_DIR}"

wget "https://sourceforge.net/projects/subread/files/subread-${VERSION}/subread-${VERSION}-Linux-x86_64.tar.gz/download" \
    -O "subread-${VERSION}-Linux-x86_64.tar.gz"

tar -xzf "subread-${VERSION}-Linux-x86_64.tar.gz"

cp -r "subread-${VERSION}-Linux-x86_64"/* "${INSTALL_DIR}/"

chmod -R a+rX "${INSTALL_DIR}"
chmod -R a+rx "${INSTALL_DIR}/bin"

echo "Subread installation complete:"
echo "${INSTALL_DIR}"
