#!/bin/bash
#SBATCH --job-name=ps_rep3mwefcomb
#SBATCH --output=logs/ps_rep3mwefcomb.out
#SBATCH --error=logs/ps_rep3mwefcomb.err
#SBATCH --partition=comp72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=14:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 3
repetition=3
traits=("1715_wave_1691" "728_wave_1071" "2335_wave_1446")

# Run the R script with the specified traits and repetition
Rscript scripts/ps_combinedrep_mwef.R --trait1 ${traits[0]} --trait2 ${traits[1]} --trait3 ${traits[2]} --repetition $repetition
