#!/bin/bash
#SBATCH --job-name=combine_all coh2 
#SBATCH --partition=tres72
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --qos comp
#SBATCH --tasks-per-node=32
#SBATCH --time=1:00:00
#SBATCH --output=logs/combine_coh2_%j.out
#SBATCH --error=logs/combine_coh2_%j.err

set -euo pipefail

mkdir -p logs

cd /home/ashmitau/GenCorr/

module load python/3.9.15  

python scripts/combine.py
