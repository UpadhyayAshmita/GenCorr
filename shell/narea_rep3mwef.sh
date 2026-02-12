#!/bin/bash
#SBATCH --job-name=narea_rep3mwef
#SBATCH --output=logs/narea_rep3mwef.out
#SBATCH --error=logs/narea_rep3mwef.err
#SBATCH --partition=comp06
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6  # Total of 6 cores (3 traits x 2 CV schemes)
#SBATCH --time=06:00:00

# Load required modules
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Define traits and repetition
repetition=3
traits=("1703_wave_1704" "2235_wave_2237" "741_wave_735")

# Run the R script with the specified traits and repetition
Rscript scripts/narea_combinedrep_mwef.R \
  --trait1 "${traits[0]}" \
  --trait2 "${traits[1]}" \
  --trait3 "${traits[2]}" \
  --repetition "$repetition"
