#!/bin/bash
#SBATCH --job-name=coh2_narea
#SBATCH --partition=comp72
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --output=logs/coh2_%j.out
#SBATCH --error=logs/coh2_%j.err

set -euo pipefail
module purge
module load R

# deterministic math (prevents BLAS oversubscription + tiny nondeterminism)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p logs

for i in $(seq 0 49); do
  start_=$((i * 43 + 350))
  end_=$(((i + 1) * 43 + 350))

  srun --exclusive -N1 -n1 Rscript ./scripts/b_coh2_scripts/coh2_block.R narea "$start_" "$end_" &

  if (( $(jobs -r | wc -l) >= 32 )); then
    wait -n
  fi
done

wait
