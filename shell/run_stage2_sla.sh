#!/bin/bash
#SBATCH -J Sla_complete_efmw
#SBATCH -p comp06
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --time=6:00:00
#SBATCH --array=1-6
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
# keep this minimal like your cloud72 jobs
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

# --- map array id -> trait_index + cv_scheme ---
# traits you mentioned earlier:
traits=("1683_wave_1666" "1640_wave_1655" "738_wave_1111")
schemes=("CV1" "CV2")

k=$((SLURM_ARRAY_TASK_ID - 1))
trait_index=$(( k / 2 + 1 ))          # 1..3
scheme_index=$(( k % 2 ))             # 0..1

trait=${traits[$((trait_index - 1))]}
cv_scheme=${schemes[$scheme_index]}

echo "Running on $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}  TASK=${SLURM_ARRAY_TASK_ID}"
echo "trait_index=${trait_index} trait=${trait} cv_scheme=${cv_scheme}"
module list

# run (edit the path to your R script file)
Rscript scripts/completemodel_parallel/rerun_sla_completemodel_mwef.R \
  --trait "${trait}" \
  --trait_index "${trait_index}" \
  --cv_scheme "${cv_scheme}"
