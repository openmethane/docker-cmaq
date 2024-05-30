#!/usr/bin/env bash

# Setup the build environment for iotools
ARCH=$(uname -m)
BLD_OS=$(/bin/uname -s)$(/bin/uname -r | cut -d. -f1)
export ARCH
export BLD_OS
export BIN=${BLD_OS}_${ARCH}gfort
export BASEDIR=$PWD
export CPLMODE=nocpl
export LD_LIBRARY_PATH="/opt/wrf/libs/lib:${LD_LIBRARY_PATH}"
export PATH="/opt/wrf/libs/bin:${PATH}"