#!/bin/bash
#
#PBS -N Run_mint_wrapper
#PBS -l select=1:ncpus=1:mem=128gb
#PBS -l walltime=48:00:00

source ~/.bashrc

cd $PBS_O_WORKDIR

conda activate Rbase-scNMT

Rscript --no-save --no-restore ./run_mint_wrapper.R
