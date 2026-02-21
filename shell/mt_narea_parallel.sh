#!/bin/bash
#SBATCH -J narea_mt_cv1
#SBATCH -p cloud72
#SBATCH -q cloud
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=09:00:00        # Slightly longer for CV1 compatibility
#SBATCH --array=1-20           # 20 reps per CV scheme
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

# Load modules
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

# --- CONFIGURATION ---
# Set this variable to "CV1" or "CV2" before submitting
# Or pass it as an argument: sbatch mt_narea_all_cvs.sh CV1
SCHEME=${1:-"CV1"}

TRAIT="narea"
WAVE1="wave_2262_wave_2275"
WAVE2="wave_846_wave_749"

echo "------------------------------------------------"
echo "Running TRAIT: $TRAIT"
echo "Running SCHEME: $SCHEME"
echo "Running REPETITION: $SLURM_ARRAY_TASK_ID"
echo "------------------------------------------------"
#
# Rscript scripts/13_mt_alltrait_parallel.R \
#     --trait "narea" \
#     --cv_scheme "$SCHEME" \
#     --w1 "$WAVE1" \
#     --w2 "$WAVE2" \
#     --rep "$SLURM_ARRAY_TASK_ID"


Rscript scripts/13_mt_alltrait_MWEF_parallel.R \
    --trait "narea" \
    --cv_scheme "$SCHEME" \
    --w1 "$WAVE1" \
    --w2 "$WAVE2" \
    --rep "$SLURM_ARRAY_TASK_ID"
