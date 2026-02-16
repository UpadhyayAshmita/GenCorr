#!/bin/bash
#SBATCH -J mt_narea_array
#SBATCH -p tres72
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32                  # Restored to 32 to satisfy QOSMinCpuNotSatisfied
#SBATCH --time=36:00:00
#SBATCH --array=0-1            # 0=CV1, 1=CV2
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

# Load modules
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2

# Map Array ID to Scheme
SCHEMES=("CV1" "CV2")
MY_SCHEME=${SCHEMES[$SLURM_ARRAY_TASK_ID]}

# Define Wavelength Columns
WAVE1="wave_2262_wave_2275"
WAVE2="wave_846_wave_749"

echo "Running Scheme: $MY_SCHEME on $SLURM_CPUS_PER_TASK cores"

Rscript scripts/mt_narea_3trait.R \
    --cv_scheme "$MY_SCHEME" \
    --w1 "$WAVE1" \
    --w2 "$WAVE2"
# Rscript scripts/13_mt_narea_3trait_mwef.R \
#     --cv_scheme "$MY_SCHEME" \
#     --w1 "$WAVE1" \
#     --w2 "$WAVE2"
