#!/bin/bash

#PBS -P lp86
#PBS -q normal
#PBS -N fwd_run
#PBS -l walltime=12:00:00,mem=96GB
#PBS -l storage=scratch/lp86+gdata/hh5
#PBS -l ncpus=16
#PBS -l wd
###PBS -j oe

source /home/563/ns0890/programs/cmaq_adj/scripts/module_cmaq.sh

###cmaq_adj crashes when environment variable set with no value.
unset PBS_NCI_IMAGE

./run.adj.fwd_d03.bnmk &> fwd.bnmk.log
