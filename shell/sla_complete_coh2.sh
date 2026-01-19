#!/bin/bash
#SBATCH --job-name=sla_coh2_complete
#SBATCH --partition=condo
#SBATCH --constraint=samuelbf
#SBATCH --account=samuelbf
#SBATCH --qos=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=72:00:00
#SBATCH --output=logs/sla_coh2_%j.out
#SBATCH --error=logs/sla_coh2_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs output

IMG="$SLURM_SUBMIT_DIR/r-ver_4.2.2.sif"
RLIB="$HOME/Rlibs_422"

# throttle: run at most N jobs at once (N=cpus requested)
max_jobs=${SLURM_CPUS_PER_TASK}

for i in {0..49}; do
  start_=$(( i*43 + 350 ))
  end_=$(( (i+1)*43 + 350 ))

  apptainer exec -e --env R_LIBS_USER="$RLIB" "$IMG" \
    Rscript scripts/coh2_scripts/coh2_block.R sla "$start_" "$end_" &

  while [ "$(jobs -r | wc -l)" -ge "$max_jobs" ]; do
    sleep 2
  done
done

wait
echo "DONE: sla coh2 blocks finished"
