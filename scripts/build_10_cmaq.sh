#!/usr/bin/env bash
# Builds CMAQ

set -x
set -e

source $PWD/scripts/common.sh

ROOT=$PWD

# Build CMAQ
cmaq_dirname=CMAQv5.0.2_notpollen
pushd ${cmaq_dirname} || exit

# Link in the required libraries
mkdir -p /opt/cmaq/${cmaq_dirname}/lib/${ARCH}/gcc/
ln -s /opt/cmaq/ioapi-3.1 /opt/cmaq/${cmaq_dirname}/lib/ioapi-3.1

pushd scripts

#  Build the builder first
pushd build
  echo "Building build"
  ./bldit.bldmake
popd

# stenex has a different named run script
pushd stenex
  echo "Building stenex"
  ./bldit.se >&! bldit.stenex.log
popd


# Loop over the other models to be built
for item in pario icon cctm; do
  pushd $item
  echo "Building $item"
  ./bldit.$item
  popd
done

pushd bcon
  echo "Building bcon"
  # Build in the profile and m3conc modes
  MODTYPE=m3conc ./bldit.bcon
  MODTYPE=profile ./bldit.bcon
popd

pushd mcip/src

# MCIP doesn't have an associated blidit script
echo "Building mcip"

FC=mpif90
FFLAGS="-g -O0 -I${ROOT}/ioapi-3.1/${BIN} $(nf-config --fflags)"
LIBS="-L${ROOT}/ioapi-3.1/${BIN} -lioapi -fopenmp $(nf-config --flibs)"

FC=$FC FFLAGS=$FFLAGS LIBS=$LIBS make

[[ -f mcip.exe ]] || { echo "MCIP failed to build"; exit 1; }
popd

# Check that BCON and ICON were built
BUILD_DIR=BLD_CH4only


ls bcon/$BUILD_DIR
ls icon/$BUILD_DIR
ls cctm/$BUILD_DIR
[[ -f bcon/$BUILD_DIR/BCON_CH4only_${BIN}_profile_CH4only ]] || { echo "BCON failed to build"; exit 1; }
[[ -f bcon/$BUILD_DIR/BCON_CH4only_${BIN}_m3conc_CH4only ]] || { echo "BCON failed to build"; exit 1; }
[[ -f icon/$BUILD_DIR/ICON_CH4only_${BIN}_profile_CH4only ]] || { echo "ICON failed to build"; exit 1; }
[[ -f cctm/$BUILD_DIR/CCTM_CH4only_${BIN} ]] || { echo "CCTM failed to build"; exit 1; }
