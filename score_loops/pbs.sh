#!/bin/bash
#PBS -l select=1:ncpus=1:mem=64G
#PBS -l walltime=01:00:00
#PBS -N get_scores
#PBS -j oe
######PBS -W depend=afterok:254600.cherry

cd $PBS_O_WORKDIR

module load conda
conda activate 3dgenome 

./get_scores_mcool.py
