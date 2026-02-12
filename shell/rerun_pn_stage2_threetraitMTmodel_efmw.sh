#!/bin/bash
#SBATCH -J plsrNarea_3synMT_efmw
#SBATCH -p comp06
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --time=6:00:00
#SBATCH --array=1-2
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

# ---- three synthetic wave ratio columns for plsr_narea (same ones you used before) ----
waves=("wave_1253_wave_376" "wave_867_wave_732" "wave_1339_wave_2266")
waves_csv=$(IFS=,; echo "${waves[*]}")   # "wave_1253_wave_376,wave_867_wave_732,wave_1339_wave_2266"

# ---- CV schemes (2 jobs) ----
schemes=("CV1" "CV2")
cv_scheme=${schemes[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Running on $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}  TASK=${SLURM_ARRAY_TASK_ID}"
echo "waves=${waves_csv}"
echo "cv_scheme=${cv_scheme}"
module list

# ---- run your EFMW 3-synthetic-trait MT model for plsr_narea ----
Rscript scripts/completemodel_parallel/stage2_pn_three_syntraitMT_model_efmw.R \
  --waves "${waves_csv}" \
  --trait_index 1 \
  --cv_scheme "${cv_scheme}"
