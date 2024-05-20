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

# List of models to build
# The first item is the directory to build
# The second item is the target
OPTIONS=build,bldmake pario,pario stenex,se icon,icon bcon,bcon mcip,mcip cctm,cctm


OLDIFS=$IFS
IFS=',';
for item in $OPTIONS; do
    # split the string using a comma into $1 and $2
    set -- $item;
    pushd $1
    echo "Building $item"
    ./bldit.$item >&! bldit.$item.log
    popd
done
IFS=$OLDIFS

popd
popd

rm ${cmaq_dirname}.zip