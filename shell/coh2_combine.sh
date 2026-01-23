#!/bin/bash
#SBATCH --job-name=coh2_combine
#SBATCH --output=logs/coh2_combine.out
#SBATCH --error=logs/coh2_combine.err
#SBATCH --partition=comp06
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6 
#SBATCH --time=01:00:00

set -euo pipefail

mkdir -p logs

cd /home/ashmitau/GenCorr/

module load python/3.9.15  

python scripts/combine.py
