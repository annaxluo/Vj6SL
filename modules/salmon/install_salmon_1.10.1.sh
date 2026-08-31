#!/usr/bin/env bash
set -euo pipefail

# Install Salmon 1.10.1
# Building on Rocky 0 requires Boost 1.82, linked statically. 

VERSION="1.10.1"
BOOST_VERSION="1.82.0"
BOOST_UNDERSCORE="1_82_0"
TOOL="salmon"

SOFTWARE_ROOT="${SOFTWARE_ROOT}"
INSTALL_DIR="${SOFTWARE_ROOT}/${TOOL}/${VERSION}"
BUILD_DIR="${TMPDIR:-/tmp}/${TOOL}-${VERSION}-build-${USER:-user}-$$"

mkdir -p "${INSTALL_DIR}"
mkdir -p "${BUILD_DIR}"

cd "${BUILD_DIR}"

wget -qO- "https://github.com/COMBINE-lab/salmon/archive/refs/tags/v${VERSION}.tar.gz" \
    | tar -xzf -

wget -qO- "https://boostorg.jfrog.io/artifactory/main/release/${BOOST_VERSION}/source/boost_${BOOST_UNDERSCORE}.tar.bz2" \
    | tar -xjf -

unset CPPFLAGS || true
unset LDFLAGS || true

echo "Building Boost ${BOOST_VERSION}"
cd "boost_${BOOST_UNDERSCORE}"
./bootstrap.sh --prefix="${INSTALL_DIR}"
./b2 cxxflags=-fPIC cflags=-fPIC link=static -a
./b2 install

cd "${BUILD_DIR}"

echo "Building Salmon ${VERSION}"
cd "salmon-${VERSION}"
mkdir -p build
cd build

cmake \
    -DNO_IPO=TRUE \
    -DFETCH_STADEN=TRUE \
    -DBOOST_ROOT="${INSTALL_DIR}" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
    ..

make -j "${THREADS:-4}"
make install

chmod -R a+rX "${INSTALL_DIR}"

echo "Salmon installation complete:"
echo "${INSTALL_DIR}"

