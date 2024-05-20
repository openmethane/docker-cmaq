#!/usr/bin/env bash
# Builds the CMAQ adjunct

set -x
set -e

cd /opt/cmaq/cmaq_adj/BLDMAKE_git
make

cd /opt/cmaq/cmaq_adj/scripts
csh bldit.adjoint.bwd.CH4only
# csh bldit.adjoint.fwd.CH4only

# TODO Verify the output exists