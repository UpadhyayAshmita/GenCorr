#!/bin/bash
#SBATCH --job-name=narea_completecomb
#SBATCH --output=logs/narea_completecomb_%A_%a.out
#SBATCH --error=logs/narea_completecomb_%A_%a.err
#SBATCH --partition=comp06
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3  # Use 3 cores per task for mclapply
#SBATCH --time=06:00:00
#SBATCH --array=1-6  # Array for 3 traits x 2 CV schemes

# Define traits and CV schemes
declare -a traits=("2262_wave_2275" "846_wave_749" "1147_wave_393")
declare -a cv_schemes=("CV2" "CV1")

# Calculate the trait and CV scheme based on SLURM_ARRAY_TASK_ID
trait_index=$(( (SLURM_ARRAY_TASK_ID + 1) / 2 ))  # 1, 1, 2, 2, 3, 3
cv_index=$(( (SLURM_ARRAY_TASK_ID - 1) % 2 + 1 ))  # 1, 2, 1, 2, 1, 2

# Select the trait and CV scheme
trait=${traits[$trait_index-1]}
cv_scheme=${cv_schemes[$cv_index-1]}

# Run the R script with the specified trait, trait index, and CV scheme
Rscript scripts/narea_completemodel.R --trait "$trait" --trait_index "$trait_index" --cv_scheme "$cv_scheme"
