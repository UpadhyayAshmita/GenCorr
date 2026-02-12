#!/bin/bash
#SBATCH --job-name=ps_rep4mwefcomb
#SBATCH --output=logs/ps_rep4mwefcomb.out
#SBATCH --error=logs/ps_rep4mwefcomb.err
#SBATCH --partition=comp72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=14:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits for repetition 4
repetition=4
traits=("1208_wave_1276" "1424_wave_1420" "2233_wave_2221")

# Run the R script with the specified traits and repetition
Rscript scripts/ps_combinedrep_mwef.R --trait1 ${traits[0]} --trait2 ${traits[1]} --trait3 ${traits[2]} --repetition $repetition
