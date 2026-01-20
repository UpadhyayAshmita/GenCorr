#!/bin/bash
#SBATCH --job-name=ps_coh2_complete 
#SBATCH --partition=comp72
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --qos comp
#SBATCH --tasks-per-node=32
#SBATCH --mem=96GB
#SBATCH --time=72:00:00
#SBATCH --output=logs/ps_coh2_%j.out
#SBATCH --error=logs/ps_coh2_%j.err



## configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

for i in {0..49}; do
 start_=$((($i * 43) + 350))
 end_=$((($i + 1) * 43 + 350))
 Rscript scripts/coh2_scripts/coh2_block.R ps $start_ $end_ &
done