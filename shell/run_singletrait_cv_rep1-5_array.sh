#!/bin/bash
#SBATCH -J pn_st_cv
#SBATCH -p comp06
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --time=06:00:00
#SBATCH --array=1-5
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge

# Load the full stack (this is the key)
module load gcc/9.3.1
module load mkl/19.0.5
module load R/4.2.2

# sanity print (goes into each .out)
echo "Node: $(hostname)"
module list
ldd "$(which R)" | egrep -i "blas|mkl|openblas|lapack" || true

REP=${SLURM_ARRAY_TASK_ID}
Rscript ./scripts/i_run_pn_singletraits_cv_rep1-5.R "${REP}"
