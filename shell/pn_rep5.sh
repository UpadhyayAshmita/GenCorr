#!/bin/bash
#SBATCH --job-name=pn_rep5
#SBATCH --output=logs/pn_rep5.out
#SBATCH --error=logs/pn_rep5.err
#SBATCH --partition=comp06
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=06:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Run the R script with arguments
Rscript scripts/pn_rep5.R --trait1 "1299_wave_738" --trait2 "1375_wave_375" --trait3 "1138_wave_744"
