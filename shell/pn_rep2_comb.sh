#!/bin/bash
#SBATCH --job-name=pn_rep2comb
#SBATCH --output=logs/pn_rep2comb.out
#SBATCH --error=logs/pn_rep2comb.err
#SBATCH --partition=comp72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=14:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 2
repetition=2
traits=("1332_wave_2422" "1698_wave_1707" "1382_wave_1791")

# Run the R script with the specified traits and repetition
Rscript scripts/pn_combinedrep.R --trait1 ${traits[0]} --trait2 ${traits[1]} --trait3 ${traits[2]} --repetition $repetition
