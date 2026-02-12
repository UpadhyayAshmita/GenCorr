#!/bin/bash
#SBATCH -J sla_W1_cv
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
module load gcc/9.3.1 mkl/19.0.5 R/4.3.0

# ---- W1 only ----
trait_index=1
trait="1683_wave_1666"

# map array id -> CV1/CV2
schemes=("CV1" "CV2")
cv_scheme=${schemes[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Running on $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}  TASK=${SLURM_ARRAY_TASK_ID}"
echo "trait_index=${trait_index} trait=${trait} cv_scheme=${cv_scheme}"
module list

Rscript scripts/completemodel_parallel/rerun_sla_completemodel_efmw.R \
  --trait "${trait}" \
  --trait_index "${trait_index}" \
  --cv_scheme "${cv_scheme}"
