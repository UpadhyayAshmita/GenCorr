#!/bin/bash
#SBATCH -J Sla_2synMT_repwise
#SBATCH -p comp06
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --array=1-20
#SBATCH --mem=0
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

# prevent oversubscription/memory spikes
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

waves=("wave_1683_wave_1666" "wave_1640_wave_1655" "wave_738_wave_1111")
waves_csv=$(IFS=,; echo "${waves[*]}")

REP="${SLURM_ARRAY_TASK_ID}"

# set scheme here, submit twice (CV1 then CV2)
CV_SCHEME="${CV_SCHEME:-CV1}"

RSCRIPT="scripts/completemodel_parallel/stage2_sla_two_syntraitMT_model_efmw.R"

echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}  Rep: ${REP}  Scheme: ${CV_SCHEME}"
echo "waves_csv='${waves_csv}'"
module list

Rscript "${RSCRIPT}" \
  --waves "${waves_csv}" \
  --trait_index 1 \
  --cv_scheme "${CV_SCHEME}" \
  --rep "${REP}"
