#!/bin/bash
#SBATCH --job-name=pnTrep4comb
#SBATCH --output=logs/pn_STrep4comb.out
#SBATCH --error=logs/pn_STrep4comb.err
#SBATCH --partition=cloud72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 4
repetition=4
trait1="667_wave_573"

# Run the R script with the specified trait and repetition
Rscript scripts/pnT_combined.R --trait1 $trait1 --repetition $repetition
