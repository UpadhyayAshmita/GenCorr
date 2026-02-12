#!/bin/bash
#SBATCH -J Narea_MT_EFMW_2syn
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

# ---- choose EXACTLY TWO synthetic wave columns ----
# replace these two with the two you want to run
waves=("wave_1511_wave_2314" "wave_1039_wave_405")
waves_csv=$(IFS=,; echo "${waves[*]}")

# ---- CV schemes (2 array tasks) ----
schemes=("CV1" "CV2")
cv_scheme="${schemes[$((SLURM_ARRAY_TASK_ID - 1))]}"

echo "Running on $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}  TASK=${SLURM_ARRAY_TASK_ID}"
echo "cv_scheme=${cv_scheme}"
echo "waves=${waves_csv}"
module list

# ---- run 2-synthetic-trait EFMW MT model ----
Rscript scripts/completemodel_parallel/stage2_narea_three_syntraitMT_model_efmw.R \
  --waves "${waves_csv}" \
  --trait_index 1 \
  --cv_scheme "${cv_scheme}"
