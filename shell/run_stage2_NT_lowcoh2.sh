#!/bin/bash
#SBATCH -J NT_lowcoh2
#SBATCH -p comp06
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --time=6:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

# ---- your single wave ratio column ----
wave_col="wave_1812_wave_1574"   # change if your column name is different for different target trait (stays same for EFMW and MWEF for same target trait)

echo "Running on $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "wave_col=${wave_col}"
module list

# This R script runs BOTH CV1 and CV2 and writes:
# NT_CV1.csv, acNT_CV1.csv, NT_CV2.csv, acNT_CV2.csv and change your path to the Rscript as you change the target trait
Rscript scripts/g_lowest_coh2/g_2_NT_mwef_lowcoh2.R \
  --wave_col "${wave_col}"
