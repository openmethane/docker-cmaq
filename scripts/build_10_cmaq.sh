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
for item in pario icon bcon mcip cctm; do
  pushd $item
  echo "Building $item"
  ./bldit.$item >&! bldit.$item.log
  popd
done

popd
popd

rm ${cmaq_dirname}.zip