#!/bin/bash
#SBATCH --job-name=ps_rep5mwefcomb
#SBATCH --output=logs/ps_rep5mwefcomb.out
#SBATCH --error=logs/ps_rep5mwefcomb.err
#SBATCH --partition=comp72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=14:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 5
repetition=5
trait3="1679_wave_1669"

# Run the R script with the specified trait and repetition
Rscript scripts/ps_combinedrep_mwef.R --trait3 $trait3 --repetition $repetition
