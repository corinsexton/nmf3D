#!/bin/bash
#PBS -l select=1:ncpus=1:mem=64G
#PBS -l walltime=50:05:00
#PBS -N create_gd
#PBS -j oe
######PBS -W depend=afterok:254600.cherry

cd $PBS_O_WORKDIR

./summarize.sh

