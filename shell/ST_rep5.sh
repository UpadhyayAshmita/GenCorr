#!/bin/bash
#SBATCH --job-name=sla_STrep5comb
#SBATCH --output=logs/sla_STrep5comb.out
#SBATCH --error=logs/sla_STrep5comb.err
#SBATCH --partition=cloud72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=15:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 4
repetition=5
trait1="915_wave_1137"

# Run the R script with the specified trait and repetition
Rscript scripts/ST_combined.R --trait1 $trait1 --repetition $repetition
