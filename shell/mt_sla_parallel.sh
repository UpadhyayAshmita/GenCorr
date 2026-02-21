#!/bin/bash
#SBATCH -J par_sla_mt_cv1
#SBATCH -p cloud72
#SBATCH -q cloud
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=09:00:00
#SBATCH --array=1-20
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

# Default to CV1 if no argument is passed
SCHEME=${1:-"CV1"}

# Wavelengths for SLA
WAVE1="wave_1226_wave_1270" #change your value based on selected synthetic trait for the target trait
WAVE2="wave_2392_wave_631"


echo "Running SLA | Scheme: $SCHEME | Rep: $SLURM_ARRAY_TASK_ID"

# Rscript scripts/13_mt_sla__parallel.R \
#     --trait "sla" \
#     --cv_scheme "$SCHEME" \
#     --w1 "$WAVE1" \
#     --w2 "$WAVE2" \
#     --rep "$SLURM_ARRAY_TASK_ID"

Rscript scripts/13_mt_sla_MWEF_parallel.R \
    --trait "sla" \
    --cv_scheme "$SCHEME" \
    --w1 "$WAVE1" \
    --w2 "$WAVE2" \
    --rep "$SLURM_ARRAY_TASK_ID"
