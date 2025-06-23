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
- From breakdown data, we can select the synthetic trait/ wave-ratios and get the phenotypic data of the target traits and the steps can be followed using the scripts in the pre-processing scripts in the scripts https://github.com/UpadhyayAshmita/GenCorr/blob/master/scripts/pre-processing.R. Also we can get the heatmaps using the scripts in the pre-processing scripts.
- After selecting three synthetic traits and getting their phenotype value, we run the script https://github.com/UpadhyayAshmita/GenCorr/blob/master/scripts/model_fit.R , in this script we first get our kinship matrix from the genotype data we have for the sorghum diversity pannel. Then we will fit the first stage model, i.e., we will get the BLUEs of the synthetic traits and target traits, respectively.
- After getting the BLUEs from the synthetic and target traits, we move forward to the second stage of the model, i.e., getting the predicted value of the genotypes. This second stage can be performed using the script below inside the folder,
https://github.com/UpadhyayAshmita/GenCorr/blob/master/scripts/second_stage.R . The second stage has three cross validation scheme i.e CV1, CV2, and single-trait for each target trait. The script includes all three scheme for four targets traits. We can obtain the accuracy and genomic estimated-breeding values for all genotypes as the output from this stage.

Here comes the end of the complete model, where all dataset is being throughout the steps.
We did the 20:80 approach as mentioned in the flowchart below:
![image](https://github.com/user-attachments/assets/bef1817e-c94c-423e-aa3b-db033f689d92)

- The sampling, i.e, random selection of the lines from the total line, was done using a script https://github.com/UpadhyayAshmita/GenCorr/blob/master/scripts/sampling_rep.R, where we can select 5 sets of phenotypes for downstream analysis. Each randomly selected line is treated as a replication. 
- Now, within each replication we selected three synthetic traits for each target trait, estimated the BLUEs for the selected synthetic traits and target traits similar to earlier complete approach and the pre-processing, selection of synthetic traits and then estimation of all 5 replication for 4 different traits can be found in the folder Replication folder with each traits R scripts.
- 
  
  

- 
