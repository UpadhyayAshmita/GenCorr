#!/bin/bash
#SBATCH --job-name=singletrait_Allrepss
#SBATCH --partition=comp72
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module purge
module load R

mkdir -p logs

# Usage:
#   sbatch shell/singletrait_secondstageallrep.sh <trait> <rep_id>
#
# Examples:
#   sbatch singletrait_allreps.sh narea 3
#   sbatch singletrait_allreps.sh sla  5
#   sbatch singletrait_allreps.sh pn   2
#   sbatch singletrait_allreps.sh ps   1

trait="${1:?Missing trait (ps|pn|narea|sla)}"
rep_id="${2:?Missing rep_id (1-5)}"

case "${trait}" in
  narea)
    rscript_path="./scripts/11_singletrait_narea_secondstageallrep.R"
    ;;
  sla)
    rscript_path="./scripts/11_singletrait_sla_secondstageallrep.R"
    ;;
  pn)
    rscript_path="./scripts/11_singletrait_pn_secondstageallrep.R"
    ;;
  ps)
    rscript_path="./scripts/11_singletrait_ps_secondstageallrep.R"
    ;;
  *)
    echo "ERROR: Unknown trait '${trait}'. Use: ps | pn | narea | sla"
    exit 1
    ;;
esac

echo "Running on $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Trait: ${trait}"
echo "Rep: ${rep_id}"
echo "Script: ${rscript_path}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

# Prevent BLAS/OpenMP oversubscription (safe default)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript "${rscript_path}" "${rep_id}"
