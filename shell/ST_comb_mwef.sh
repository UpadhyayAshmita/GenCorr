#!/bin/bash
#SBATCH --job-name=STcombinedmwef
#SBATCH --output=logs/STcombinedmwef_%A_%a.out
#SBATCH --error=logs/STcombinedmwef_%A_%a.err
#SBATCH --partition=cloud72
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --array=1-10  # Array for 5 repetitions x 2 CV schemes

# Configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

# Array of traits
declare -a traits=("367_wave_676" "618_wave_1713" "717_wave_490" "1954_wave_2292" "915_wave_1137")
declare -a cv_schemes=("CV2" "CV1")

# Calculate the repetition and CV scheme based on SLURM_ARRAY_TASK_ID
rep_index=$(( (SLURM_ARRAY_TASK_ID + 1) / 2 ))  # 1, 1, 2, 2, 3, 3, etc.
cv_index=$(( (SLURM_ARRAY_TASK_ID - 1) % 2 + 1 ))  # Alternates between 1 and 2

trait1=${traits[$rep_index-1]}
cv_scheme=${cv_schemes[$cv_index-1]}

# Run the R script with the specified trait, repetition, and CV scheme
Rscript scripts/ST_combined_mwef.R --trait1 $trait1 --repetition $rep_index --cv_scheme $cv_scheme
