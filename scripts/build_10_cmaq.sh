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
# unzip CMAQ-5.0.2.zip
# pushd CMAQ-5.0.2 || exit

