#!/bin/bash
#SBATCH --job-name=pn_rep3
#SBATCH --output=logs/pn_rep3.out
#SBATCH --error=logs/pn_rep3.err
#SBATCH --partition=comp06
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=06:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Run the R script with arguments
Rscript scripts/pn_rep3.R --trait1 "1133_wave_742" --trait2 "1703_wave_1706" --trait3 "1248_wave_1203"
