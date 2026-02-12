#!/bin/bash
#SBATCH --job-name=pnT_efmw
#SBATCH --output=logs/pnTefmw_%A_%a.out
#SBATCH --error=logs/pnTefmw_%A_%a.err
#SBATCH --constraint=samuelbf
#SBATCH --partition=condo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --array=1-10

# Load required modules
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Array of traits corresponding to each repetition
declare -a traits=("351_wave_1389" "2326_wave_1956" "1289_wave_759" "667_wave_573" "2278_wave_2111")
declare -a cv_schemes=("CV1" "CV2")

# Calculate the repetition and CV scheme index based on SLURM_ARRAY_TASK_ID
rep_index=$(( (SLURM_ARRAY_TASK_ID + 1) / 2 ))  # 1 to 5, corresponding to repetitions 1 to 5
cv_index=$(( (SLURM_ARRAY_TASK_ID - 1) % 2 ))   # 0 for CV1, 1 for CV2

# Select the correct trait and CV scheme for this task
trait=${traits[$((rep_index - 1))]}   # Use $((rep_index - 1)) to index the trait array
cv_scheme=${cv_schemes[$cv_index]}

# Run the R script with the selected trait, repetition, and CV scheme
Rscript scripts/pnT_combined.R --trait "$trait" --repetition "$rep_index" --cv_scheme "$cv_scheme"
