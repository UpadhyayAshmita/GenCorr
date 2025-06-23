#!/bin/bash

#SBATCH --job-name=SW3_rep3mwef
#SBATCH --output=logs/SW3_rep3mwefout.txt
#SBATCH --partition=tres72
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --time=12:00:00

## configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

## run
## Create multitrait prediction:
Rscript scripts/MWEF/SW3_rep3_mwef.R
