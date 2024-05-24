#!/usr/bin/env bash

# Setup the build environment for iotools
ARCH=$(uname -i)
export ARCH
export BIN=Linux2_${ARCH}gfort
export BASEDIR=$PWD
export CPLMODE=nocpl
export LD_LIBRARY_PATH="/opt/wrf/libs/lib:${LD_LIBRARY_PATH}"
export PATH="/opt/wrf/libs/bin:${PATH}"