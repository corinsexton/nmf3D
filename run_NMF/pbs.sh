#!/bin/bash
#PBS -l select=1:ncpus=4:mem=120G
#PBS -l walltime=72:00:00
#PBS -N nmf4
#PBS -j oe
######PBS -W depend=afterok:254600.cherry

cd $PBS_O_WORKDIR

time ./nmf.R

