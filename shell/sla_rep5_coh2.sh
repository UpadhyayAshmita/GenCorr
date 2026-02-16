#!/bin/bash
#SBATCH --job-name=sla_coh2_rep5
#SBATCH --partition=comp72
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --qos comp
#SBATCH --tasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH --output=logs/sla_coh2_%j.out
#SBATCH --error=logs/sla_coh2_%j.err


## configs
module purge
module load gcc/9.3.1 mkl/19.0.5 R/4.2.2 vcftools/0.1.15 plink/5.2

for i in {0..49}; do
 start_=$((($i * 43) + 350))
 end_=$((($i + 1) * 43 + 350))
 Rscript scripts/02_coh2_block_rep5.R sla $start_ $end_ &
done
