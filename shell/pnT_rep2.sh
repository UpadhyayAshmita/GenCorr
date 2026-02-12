#!/bin/bash
#SBATCH --job-name=pnTrep2comb
#SBATCH --output=logs/pn_STrep2comb.out
#SBATCH --error=logs/pn_STrep2comb.err
#SBATCH --partition=cloud72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 2
repetition=2
trait1="2326_wave_1956"

# Run the R script with the specified trait and repetition
Rscript scripts/pnT_combined.R --trait1 $trait1 --repetition $repetition
