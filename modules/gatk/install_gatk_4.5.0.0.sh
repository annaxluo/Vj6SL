#!/usr/bin/env bash
set -euo pipefail

# Install GATK 4.5.0.0 and its recommended conda environment.
#
# Usage:
#   export SOFTWARE_ROOT=/path/to/shared/software 
#   bash modules/gatk/install_gatk_4.5.0.0.sh

VERSION="4.5.0.0"
INSTALL_DIR="${SOFTWARE_ROOT}/gatk/${VERSION}"

mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

wget https://github.com/broadinstitute/gatk/releases/download/${VERSION}/gatk-${VERSION}.zip
unzip gatk-${VERSION}.zip
rm gatk-${VERSION}.zip
mv gatk-${VERSION} gatk

export PYTHONNOUSERSITE=1

# set up GATK conda environment from the YAML included with GATK.
module load conda/3-23.3.1
conda activate mamba_env
mamba env create -p "${INSTALL_DIR}/gatk_env" -f "${INSTALL_DIR}/gatk/gatkcondaenv.yml"

chmod 775 -R "${INSTALL_DIR}"
chmod 555 -R "${INSTALL_DIR}/gatk_env"

echo "GATK ${VERSION} installed in:"
echo "${INSTALL_DIR}"
echo
echo "GATK executable:"
echo "${INSTALL_DIR}/gatk/gatk"
