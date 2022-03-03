#!/bin/bash
#PBS -l select=1:ncpus=1:mem=120G
#PBS -l walltime=100:30:00
#PBS -N 4DNFI9E222YJ_peakachu_loops
#PBS -j oe
######PBS -W depend=afterok:254600.cherry

cd $PBS_O_WORKDIR

#fastq-dump SRR851838

#
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda create -n 3dgenome python=3.6 scikit-learn=0.20.2 numpy scipy pandas h5py cooler
#source activate 3dgenome
#pip install hic-straw
#git clone https://github.com/tariks/peakachu
#cd peakachu
#python setup.py install
#

sample=4DNFIJLK5WML
#sample=4DNFI9E222YJ

module load conda
conda activate 3dgenome
#peakachu depth -p raw_data/4DNFIJLK5WML.mcool::/resolutions/5000
#peakachu depth -p raw_data/4DNFI9E222YJ.mcool::/resolutions/5000

#echo "SCORE_GENOME"
time peakachu score_genome -r 5000 --balance -p raw_data/${sample}.mcool::/resolutions/5000 \
	-O call_loops/${sample}_scores_5kb.loops -m util/high-confidence.150million.5kb.pkl
echo "POOL"
#for i in call_loops/${sample}_scores_1kb/*; do peakachu pool -i $i -r 1000 -t .8 > ${i}.1kb.loops.txt; done
peakachu pool -r 5000 -i call_loops/${sample}_scores_5kb.loops -o call_loops/${sample}_scores_5kb.loops.final
