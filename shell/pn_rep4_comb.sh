#!/bin/bash
#SBATCH --job-name=pn_rep4comb
#SBATCH --output=logs/pn_rep4comb.out
#SBATCH --error=logs/pn_rep4comb.err
#SBATCH --partition=comp72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=14:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 4
repetition=4
traits=("1708_wave_1709" "1390_wave_1862" "1139_wave_1273")

# Run the R script with the specified traits and repetition
Rscript scripts/pn_combinedrep.R --trait1 ${traits[0]} --trait2 ${traits[1]} --trait3 ${traits[2]} --repetition $repetition
