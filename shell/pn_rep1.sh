#!/bin/bash
#SBATCH --job-name=pn_rep1
#SBATCH --output=logs/pn_rep1.out
#SBATCH --error=logs/pn_rep1.err
#SBATCH --partition=comp06
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=06:00:00

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Run the R script with arguments
Rscript scripts/pn_rep1.R --trait1 "753_wave_750" --trait2 "1357_wave_1800" --trait3 "758_wave_757"
