#!/usr/bin/env bash
# Builds the CMAQ adjunct

set -x
set -e

cd /opt/cmaq/cmaq_adj/BLD_fwd_CH4only
make

cd /opt/cmaq/cmaq_adj/BLD_bwd_CH4only
make
