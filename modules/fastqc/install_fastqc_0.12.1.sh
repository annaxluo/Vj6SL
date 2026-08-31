#!/usr/bin/env bash
set -euo pipefail

## Install FastQC 0.12.1 
##
## Usage:
##   export SOFTWARE_ROOT=/path/to/shared/software
##   bash modules/fastqc/install_fastqc_0.12.1.sh

VERSION="0.12.1"
TOOL="FastQC"

SOFTWARE_ROOT="${SOFTWARE_ROOT}"
INSTALL_DIR="${SOFTWARE_ROOT}/${TOOL}/${VERSION}"
BUILD_DIR="${TMPDIR:-/tmp}/fastqc-${VERSION}-build-${USER:-user}-$$"

mkdir -p "${INSTALL_DIR}"
mkdir -p "${BUILD_DIR}"

cd "${BUILD_DIR}"

wget "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${VERSION}.zip" \
    -O "fastqc_v${VERSION}.zip"

unzip "fastqc_v${VERSION}.zip"

cp -r FastQC/* "${INSTALL_DIR}/"

chmod a+rx "${INSTALL_DIR}/fastqc"
chmod -R a+rX "${INSTALL_DIR}"

echo "FastQC installed to: "
echo "${INSTALL_DIR}"
