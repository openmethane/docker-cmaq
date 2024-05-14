#!/usr/bin/env bash
# Builds CMAQ

set -x
set -e

IOAPI_VERSION="3.1"
CMAQ_VERSION=${CMAQ_VERSION:-"5.0.2"}

# Setup the build environment for iotools
export BIN=${BIN:-Linux2_x86_64gfort}
export BASEDIR=$PWD
export CPLMODE=nocpl

# Build IOAPI
wget -nv https://zenodo.org/records/1079898/files/CMAQ-5.0.2.zip  -O CMAQ-5.0.2.zip
unzip CMAQ-5.0.2.zip
pushd CMAQ-5.0.2 || exit

# Link in the required libraries
cp /opt/cmaq/templates/CMAQ/scripts/config.cmaq scripts/config.cmaq
mkdir -p /opt/cmaq/CMAQ-5.0.2/lib/${ARCH}/gcc/
ln -s /opt/cmaq/ioapi-3.1 /opt/cmaq/CMAQ-5.0.2/lib/${ARCH}/gcc/ioapi_3.1

pushd scripts

# Update the netcdf libraries
find ./ -type f -exec sed -i -e 's/-lnetcdf/-lnetcdf -lnetcdff/g' {} \;

pushd build
./bldit.bldmake >&! bldit.bldmake.log
popd

pushd stenex
./bldit.se >&! bldit.se.log
popd

pushd pario
./bldit.pario >&! bldit.pario.log
popd


popd
popd

rm CMAQ-5.0.2.zip