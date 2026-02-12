#!/bin/bash
#SBATCH --job-name=pnTrep1comb
#SBATCH --output=logs/pn_STrep1comb.out
#SBATCH --error=logs/pn_STrep1comb.err
#SBATCH --partition=cloud72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 1
repetition=1
trait1="351_wave_1389"

# Run the R script with the specified trait and repetition
Rscript scripts/pnT_combined.R --trait1 $trait1 --repetition $repetition
