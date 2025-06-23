## General

This repository contains all the scripts and data used for the GenCorr project.

## Content of Repository

-Raw data used in the study from the WEST sorghum diversity panel. To run the pipeline, we must have phenotype and genotype data, where phenotype data contains hyperspectral data and 
trait of interest data.

## Methodology
-The first step is to calculate the co-heritability, which can be done using the bash script and R script located in the co-heritability folder within the scripts.
https://github.com/UpadhyayAshmita/GenCorr/blob/master/scripts/coh2_scripts/coh2_block.R 
for all four main target traits.
- We can obtain the breakdown data for each trait using the Python scripts https://github.com/UpadhyayAshmita/GenCorr/blob/master/scripts/coh2_scripts/combine.ipynb
scripts inside the coh2 folder in the scripts.
- From breakdown data, we can select the synthetic trait/ wave-ratios and get the phenotypic data of the target traits and the steps can be followed using the scripts in the
  pre-processing scripts in the scripts scripts/pre-processing.R
  also we can get the heatmaps using the scripts in the pre-processing scripts.
- 
-  
- 
- Then, we can go for genomic prediction where we can follow two steps model and this steps can be performed using the scripts from the 

- 
