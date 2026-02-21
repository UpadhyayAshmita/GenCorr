#!/bin/bash
#SBATCH -J ps_mt_cv
#SBATCH -p cloud72
#SBATCH -q cloud
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=10:00:00
#SBATCH --array=1-20
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

# Change to CV1 or CV2 when submitting
SCHEME=${1:-"CV1"}

# Use the wavelength columns selected for PS
WAVE1="wave_2372_wave_1875"
WAVE2="wave_389_wave_1323"

echo "Running PS | Scheme: $SCHEME | Rep: $SLURM_ARRAY_TASK_ID"

# Rscript scripts/13_mt_alltrait_parallel.R \
#     --trait "plsr_sla" \
#     --cv_scheme "$SCHEME" \
#     --w1 "$WAVE1" \
#     --w2 "$WAVE2" \
#     --rep "$SLURM_ARRAY_TASK_ID"


Rscript scripts/13_mt_alltrait_MWEF_parallel.R \
    --trait "plsr_sla" \
    --cv_scheme "$SCHEME" \
    --w1 "$WAVE1" \
    --w2 "$WAVE2" \
    --rep "$SLURM_ARRAY_TASK_ID"
