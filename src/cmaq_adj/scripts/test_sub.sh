#!/bin/bash
#PBS -P q90
#PBS -q normal
#PBS -N test_run
#PBS -l walltime=20:00,mem=32GB
#PBS -l ncpus=16
#PBS -l wd

source /home/563/spt563/mods/module_cmaq.sh

#./tmp.run.fwd &> fwd.tmp.log
./tmp.run.fwd
