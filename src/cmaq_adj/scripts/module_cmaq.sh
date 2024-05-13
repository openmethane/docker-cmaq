#!/bin/bash

module purge
module load pbs
module load intel-compiler/2019.3.199
module load openmpi/4.0.2
module load hdf5/1.10.5
module load netcdf/4.7.1
module load python3/3.8.5

#NOTE: already locally installed (via pip) netCDF4 and pyproj