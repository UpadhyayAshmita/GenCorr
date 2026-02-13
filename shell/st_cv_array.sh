#!/bin/bash
#SBATCH -J ST_CV
#SBATCH -p comp06
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --time=6:00:00
#SBATCH --array=1-8
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load dnnl-cpu-iomp/2023.0.0
module load oneapi24/advisor/2023.2.0

module load R/4.3.0
module list
echo "R path: $(which Rscript)"
ldd "$(which R)" 2>/dev/null | grep -i -E "iomp|omp|mkl" || true

# Map array index -> trait + scheme
traits=(narea sla plsr_narea plsr_sla)
schemes=(EFMW MWEF)

i=$((SLURM_ARRAY_TASK_ID - 1))
trait=${traits[$(( i % 4 ))]}
scheme=${schemes[$(( i / 4 ))]}

echo "Running trait=${trait} scheme=${scheme} on $(hostname)"
Rscript scripts/e_run_singletraits_cv_completemodel.R "${trait}" "${scheme}"
