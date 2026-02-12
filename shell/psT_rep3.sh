#!/bin/bash
#SBATCH --job-name=ps_STrep3comb
#SBATCH --output=logs/ps_STrep3comb.out
#SBATCH --error=logs/ps_STrep3comb.err
#SBATCH --partition=cloud72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 3
repetition=3
trait1="467_wave_2178"

# Run the R script with the specified trait and repetition
Rscript scripts/psT_combined.R --trait1 $trait1 --repetition $repetition
