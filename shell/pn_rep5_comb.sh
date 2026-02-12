#!/bin/bash
#SBATCH --job-name=pn_rep5comb
#SBATCH --output=logs/pn_rep5comb.out
#SBATCH --error=logs/pn_rep5comb.err
#SBATCH --partition=comp72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=14:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 4
repetition=5
traits=("1299_wave_738" "1375_wave_375" "1138_wave_744")

# Run the R script with the specified traits and repetition
Rscript scripts/pn_combinedrep.R --trait1 ${traits[0]} --trait2 ${traits[1]} --trait3 ${traits[2]} --repetition $repetition
