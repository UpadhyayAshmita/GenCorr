#!/bin/bash
#SBATCH --job-name=coh2_combine
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=cloud72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00

set -euo pipefail

cd /home/ashmitau/GenCorr/
mkdir -p logs

module purge
module load python/3.9.15

which python
python -V

python scripts/combine.py
