#!/bin/bash
#SBATCH --job-name=pn_CV1_rep4
#SBATCH --output=logs/pn_CV1_rep4.out
#SBATCH --error=logs/pn_CV1_rep4.err
#SBATCH --partition=comp06
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=06:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Run the R script with arguments
Rscript scripts/pn_CV1_rep4.R --trait1 "1708_wave_1709" --trait2 "1390_wave_1862" --trait3 "1139_wave_1273"
