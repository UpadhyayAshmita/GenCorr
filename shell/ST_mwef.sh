#!/bin/bash
#SBATCH --job-name=ST_mwef
#SBATCH --output=logs/ST_mwef%A_%a.out  # %A for array job ID, %a for array task ID
#SBATCH --error=logs/ST_mwef%A_%a.err   # %A for array job ID, %a for array task ID
#SBATCH --partition=cloud72
#SBATCH --nodes=1                   # Use only 1 node
#SBATCH --ntasks=1                  # Each SLURM array task runs as a separate task
#SBATCH --cpus-per-task=1           # Each task gets 1 CPU core
#SBATCH --time=10:00:00             # Set a time limit
#SBATCH --array=1-10                # Create an array with only tasks 9 and 10 (5th repetition)

# Load required modules
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Array of traits corresponding to each repetition
declare -a traits=("367_wave_676" "618_wave_1713" "717_wave_490" "1954_wave_2292" "915_wave_1137")
declare -a cv_schemes=("CV1" "CV2")

# Calculate the repetition and CV scheme index based on SLURM_ARRAY_TASK_ID
rep_index=$(( (SLURM_ARRAY_TASK_ID + 1) / 2 ))  # 1 to 5, corresponding to repetitions 1 to 5
cv_index=$(( (SLURM_ARRAY_TASK_ID - 1) % 2 ))   # 0 for CV1, 1 for CV2

# Select the correct trait and CV scheme for this task
trait=${traits[$((rep_index - 1))]}   # Use $((rep_index - 1)) to index the trait array
cv_scheme=${cv_schemes[$cv_index]}

# Run the R script with the selected trait, repetition, and CV scheme
Rscript scripts/ST_combined_mwef.R --trait "$trait" --repetition "$rep_index" --cv_scheme "$cv_scheme"
