#!/bin/bash
#SBATCH -J ps_st_cv
#SBATCH -p comp06
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --time=6:00:00
#SBATCH --array=1-5
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

REP=${SLURM_ARRAY_TASK_ID}

Rscript scripts/i_run_ps_singletraits_cv_re1-5_trycatch.R "${REP}"
