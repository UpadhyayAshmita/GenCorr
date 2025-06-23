#!/bin/bash

#SBATCH --job-name=SW2_CV1rep1_mwef
#SBATCH --output=logs/SW2_CV1_rep1mwefout.txt
#SBATCH --partition=tres72
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --time=12:00:00

## configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

## run
## Create multitrait prediction:
Rscript scripts/MWEF/SW2_CV1_rep1_mwef.R
