#!/bin/bash
#SBATCH --job-name=narea_coh2_rep2
#SBATCH --partition=comp72
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --qos comp
#SBATCH --tasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH --output=logs/narea_coh2_%j.out
#SBATCH --error=logs/narea_coh2_%j.err




## configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

for i in {0..49}; do
 start_=$((($i * 43) + 350))
 end_=$((($i + 1) * 43 + 350))
 Rscript scripts/b_coh2_scripts/coh2_block_rep2.R narea $start_ $end_ &
done
