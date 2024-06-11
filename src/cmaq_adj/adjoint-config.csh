#!/bin/csh

# CMAQ Adjoint configuration script
# Common configuration for the fwd and bwd adjoints
# This has been tested using gcc9 and openmpi

# Code locations

setenv M3MODEL /opt/cmaq/cmaq_adj
setenv CMAQ_HOME /opt/cmaq/CMAQv5.0.2_notpollen
setenv IOAPI_HOME /opt/cmaq/ioapi-3.1

#> architecture & compiler specific settings

setenv system "`/bin/uname -m`"
setenv BLD_OS "`/bin/uname -s``/bin/uname -r | cut -d. -f1`"

# CHANGE: set M3LIB to base directory for MPI, netCDF, IOAPI, STENEX, PARIO libraries
setenv lib_basedir ${CMAQ_HOME}/lib
setenv COMPILER gcc
setenv compiler_ext gfort
setenv M3LIB ${lib_basedir}/${system}/{$COMPILER}

# Compilers
setenv FC mpif90
setenv CC mpicc
setenv FP $FC

# CHANGE: Set location of MPICH if using multiple processors
setenv MPICH  "/opt/venv"

# CHANGE: Set location for stenex library/include/and mod files
# setenv STENEX   $M3LIB/se_noop
setenv STENEX  ${M3LIB}/se_snl

# CHANGE: Set compiler flags
# These have been tested with gfortran9
# -march=native -mtune=native
setenv F_FLAGS "-ffixed-line-length-none -fd-lines-as-comments -fallow-argument-mismatch -dI -cpp -I. -g -O2"
setenv CPP_FLAGS ""
setenv C_FLAGS    "-g -DFLDMN -I${MPICH}/include"
setenv LINK_FLAGS "-fopenmp"
# CHANGE: Set location of libraries/include files
setenv IOAPI_FLAGS  "-L${IOAPI_HOME}/${BLD_OS}_${system}${compiler_ext} -lioapi"
setenv NETCDF_FLAGS  "`nf-config --flibs`"
setenv PARIO_FLAGS "-L${M3LIB}/pario -lpario"
setenv MPICH_FLAGS "-L${MPICH}/lib -lmpich"

# Set location of M3Bld executable
setenv Blder  "$M3MODEL/BLDMAKE_git/bldmake -verbose"

setenv ICL_IOAPI  "${IOAPI_HOME}/ioapi/fixed_src"
