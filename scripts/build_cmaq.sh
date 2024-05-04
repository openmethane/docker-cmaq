#!/usr/bin/env bash
# Builds ioapi and CMAQ using the prebuilt dependencies from the WRF container

set -x
set -e

IOAPI_VERSION="3.2"
CMAQ_VERSION=${-CMAQ_VERSION:"5.0.2"}
DIR=/opt/wrf/libs

# Link to the compiled dependencies from the WRF container
export PATH=$DIR/bin:$PATH
export LDFLAGS=-L$DIR/lib
export CPPFLAGS=-I$DIR/include

# Build IOAPI
wget -nv https://www.cmascenter.org/ioapi/download/ioapi-${IOAPI_VERSION}.tar.gz  -O ioapi-v${IOAPI_VERSION}.tar.gz
mkdir ioapi-v${IOAPI_VERSION}
tar -xzvf  ioapi-v${IOAPI_VERSION}.tar.gz -C ioapi-v${IOAPI_VERSION}
pushd ioapi-v${IOAPI_VERSION} || exit

# Copy in Makefile template
# This uses the correct configuration for the built NetCDF library
cp /opt/wrf/build/templates/ioapi/Makefile Makefile
cp ioapi/Makefile.nocpl ioapi/Makefile # This is overridden by the .nocpl.sed template
cp m3tools/Makefile.nocpl m3tools/Makefile  # This is overridden by the .nocpl.sed template

# Setup the build environment
export BIN=Linux2_x86_64gfort10
export BASEDIR=$PWD
export CPLMODE=nocpl

make

# build CMAQ
echo "Building CMAQ v${CMAQ_VERSION}"