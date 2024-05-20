#!/usr/bin/env bash
# Builds CMAQ

set -x
set -e

source $PWD/scripts/common.sh

CMAQ_VERSION=${CMAQ_VERSION:-"5.0.2"}
DOI="1079898" # This is implicitly the DOI for CMAQ 5.0.2

# Setup the build environment for iotools
ARCH=$(uname -i)
export ARCH
export BIN=Linux2_${ARCH}gfort
export BASEDIR=$PWD
export CPLMODE=nocpl

ROOT=$PWD

# Build IOAPI
cmaq_dirname=CMAQ-${CMAQ_VERSION}
wget -nv https://zenodo.org/records/${DOI}/files/${cmaq_dirname}.zip  -O ${cmaq_dirname}.zip
unzip ${cmaq_dirname}.zip
pushd ${cmaq_dirname} || exit

# Link in the required libraries
cp /opt/cmaq/templates/CMAQ/scripts/config.cmaq scripts/config.cmaq
mkdir -p /opt/cmaq/${cmaq_dirname}/lib/${ARCH}/gcc/
ln -s /opt/cmaq/ioapi-3.1 /opt/cmaq/${cmaq_dirname}/lib/${ARCH}/gcc/ioapi_3.1

pushd scripts

# Update the netcdf libraries
find ./ -type f -exec sed -i -e 's/-lnetcdf/-lnetcdf -lnetcdff/g' {} \;

#  Build the builder first
pushd build
  echo "Building build"
  ./bldit.bldmake >&! bldit.build.log
popd

# stenex has a different named run script
pushd stenex
  echo "Building stenex"
  ./bldit.se >&! bldit.stenex.log
popd


# Loop over the other models to be built
for item in pario icon bcon; do
  pushd $item
  echo "Building $item"
  ./bldit.$item
  popd
done

pushd mcip/src

# MCIP doesn't have an associated blidit script
echo "Building mcip"

FC=mpif90
FFLAGS="-O3 -I${ROOT}/ioapi-3.1/${BIN} $(nf-config --fflags)"
LIBS="-L${ROOT}/ioapi-3.1/${BIN} -lioapi -lnetcdf -lnetcdff -fopenmp"

FC=$FC FFLAGS=$FFLAGS LIBS=$LIBS make

[[ -f mcip.exe ]] || { echo "MCIP failed to build"; exit 1; }

rm $ROOT/${cmaq_dirnasme}.zip