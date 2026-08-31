#!/usr/bin/env bash
set -euo pipefail

# Install Trimmomatic 0.39.
#
# Usage:
#   export SOFTWARE_ROOT=/path/to/shared/software
#   bash modules/trimmomatic/install_trimmomatic_0.39.sh

VERSION="0.39"
INSTALL_DIR="${SOFTWARE_ROOT}/Trimmomatic/${VERSION}"
BUILD_DIR="${INSTALL_DIR}/build"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${VERSION}.zip
unzip Trimmomatic-${VERSION}.zip

mkdir -p "${INSTALL_DIR}"
cp -r Trimmomatic-${VERSION}/* "${INSTALL_DIR}/"

chmod -R a+rX "${INSTALL_DIR}"

echo "Trimmomatic ${VERSION} installed in:"
echo "${INSTALL_DIR}"
echo "JAR file:"
echo "${INSTALL_DIR}/trimmomatic-${VERSION}.jar"
