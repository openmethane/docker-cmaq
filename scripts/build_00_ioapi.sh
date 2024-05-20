#!/usr/bin/env bash
# Builds ioapi

set -x
set -e

source $PWD/scripts/common.sh

IOAPI_VERSION="3.1"


# Build IOAPI
wget -nv https://www.cmascenter.org/ioapi/download/ioapi-${IOAPI_VERSION}.tar.gz  -O ioapi-${IOAPI_VERSION}.tar.gz
mkdir ioapi-${IOAPI_VERSION}
tar -xzvf  ioapi-${IOAPI_VERSION}.tar.gz -C ioapi-${IOAPI_VERSION}
pushd ioapi-${IOAPI_VERSION} || exit

# Copy in Makefile template
# This uses the correct configuration for the built NetCDF library
cp /opt/cmaq/templates/ioapi/Makefile Makefile
cp /opt/cmaq/templates/ioapi/Makeinclude.* ioapi/
# The following makefiles will be overridden by the .nocpl.sed templates
# but need to exist for the configure step to work
#touch ioapi/Makefile m3tools/Makefile

mkdir $BIN

make configure
make

popd

# Clean up
rm ioapi-${IOAPI_VERSION}.tar.gz
