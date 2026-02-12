#!/bin/bash
#SBATCH --job-name=pnTrep3comb
#SBATCH --output=logs/pn_STrep3comb.out
#SBATCH --error=logs/pn_STrep3comb.err
#SBATCH --partition=cloud72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 3
repetition=3
trait1="1289_wave_759"

# Run the R script with the specified trait and repetition
Rscript scripts/pnT_combined.R --trait1 $trait1 --repetition $repetition
