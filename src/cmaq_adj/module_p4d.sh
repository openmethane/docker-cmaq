#!/bin/bash

#cmaq-stuff
module purge
module load pbs
module load intel-compiler/2019.3.199
module load openmpi/4.0.3
module load netcdf/4.7.1
module load hdf5/1.10.5

#python-stuff
module load python3/3.7.4

#NOTE: already locally installed (via pip) netCDF4 and pyproj
