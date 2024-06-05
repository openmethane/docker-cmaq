#!/usr/bin/env bash
# Builds the CMAQ adjunct

set -x
set -e

source $PWD/scripts/common.sh

cd /opt/cmaq/cmaq_adj/BLDMAKE_git
make

cd /opt/cmaq/cmaq_adj/scripts

# Build the Backward adjoint model
csh bldit.adjoint.bwd.CH4only

if [[ ! -f /opt/cmaq/cmaq_adj/BLD_bwd_CH4only/ADJOINT_BWD ]]; then
    echo "Backward adjoint failed to build"
    exit 1
fi

# Build the Forward adjoint model
csh bldit.adjoint.fwd.CH4only

if [[ ! -f /opt/cmaq/cmaq_adj/BLD_fwd_CH4only/ADJOINT_FWD ]]; then
    echo "Forward adjoint failed to build"
    exit 1
fi

/opt/cmaq/cmaq_adj/BLD_fwd_CH4only/ADJOINT_FWD 2>&1 | tee -a adj_fwd.log