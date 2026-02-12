#!/bin/bash
#SBATCH -J plsrSLA_3synMT_efmw
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

# ---- three synthetic wave ratio columns for plsr_sla (same ones you used before) ----
waves=("wave_2242_wave_1540" "wave_2156_wave_1385" "wave_393_wave_779")
waves_csv=$(IFS=,; echo "${waves[*]}")   # "wave_2242_wave_1540,wave_2156_wave_1385,wave_393_wave_779"

# ---- CV schemes (2 jobs) ----
schemes=("CV1" "CV2")
cv_scheme=${schemes[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Running on $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}  TASK=${SLURM_ARRAY_TASK_ID}"
echo "waves=${waves_csv}"
echo "cv_scheme=${cv_scheme}"
module list

# ---- run your EFMW 3-synthetic-trait MT model for plsr_sla ----
Rscript scripts/completemodel_parallel/rerun_ps_stage2_three_syntraitMT_model_efmw.R \
  --waves "${waves_csv}" \
  --trait_index 1 \
  --cv_scheme "${cv_scheme}"
