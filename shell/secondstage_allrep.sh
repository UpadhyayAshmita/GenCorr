#!/bin/bash
#SBATCH --job-name=mt_submit_repmodel
#SBATCH --partition=comp72
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=logs/mt_%j.out
#SBATCH --error=logs/mt_%j.err

set -euo pipefail

module purge
module load R

mkdir -p logs

# Usage:
# sbatch mt_submit.sbatch <target> <trait1> <trait2> <trait3> <rep>
target="${1:?Missing target (ps|pn|narea|sla)}"
trait1="${2:?Missing trait1}"
trait2="${3:?Missing trait2}"
trait3="${4:?Missing trait3}"
rep="${5:?Missing repetition}"

case "${target}" in
  ps)
    rscript_path="./scripts/12_ps_combinedrep.R"
    ;;
  pn)
    rscript_path="./scripts/12_pn_combinedrep.R"
    ;;
  narea)
    rscript_path="./scripts/12_narea_combinedrep.R"
    ;;
  sla)
    rscript_path="./scripts/12_sla_combinedrep.R"
    ;;
  *)
    echo "ERROR: Unknown target '${target}'. Use: ps | pn | narea | sla"
    exit 1
    ;;
esac

echo "Running on $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Target: ${target}"
echo "Traits: ${trait1}, ${trait2}, ${trait3}"
echo "Rep: ${rep}"
echo "Script: ${rscript_path}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

# Prevent BLAS/OpenMP oversubscription (important with mclapply)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript "${rscript_path}" \
  --trait1 "${trait1}" \
  --trait2 "${trait2}" \
  --trait3 "${trait3}" \
  --repetition "${rep}"
