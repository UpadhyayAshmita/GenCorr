## How to reproduce the results
- High-performance computing with high memory and enough to submit the job as per the memory and resource allocated in a bash script in shell dir
- Also, some interactive R session is necessary to pre-process the data
## Clone the dir and create this dir structure 
- After cloning the dir, phenotypic data, which is available in _ needs to be downloaded and put it inside data dir along with the genmoic data, which is also available here _, and Names_WEST_SF.csv data for filtering later on the steps needs to be downloaded.
- Dir structure should be like this:
```
GenCorr/
├── data/                    # phenotypic + genomic data + Names_WEST.csv file 
├── figure/                  #for result figure 
├── function/                #custom functions
├── logs/                    #for checking job logs
├── output/                  #for output of pre-processing, first stage and second stage model
├── output_complete/         #for coh2_block complete model output 
├── output1/                 #coh2_block replication 1 for subset 1 output
├── output2/                 #coh2_block replication 2 for subset 2 output
├── output3/                 #coh2_block replication 3 for subset 3 output
├── output4/                 #coh2_block replication 4 for subset 4 output
├── output5/                 #coh2_block replication 5 for subset 5 output
├── scripts/                 # for R scripts 
└── shell/                   # for SLURM scheduler or bash script 
```



## For partitioning the complete vs replication 1-5 dataset before calculating coheritability 
Run through an interactive R session 
```
srun --partition=batch --ntasks=1 --cpus-per-task=32 --mem=164G --time=5:00:00 --pty bash
```
```
./scripts/a_sampling_rep.R 
```
## Coheritability ratios generation for the complete dataset and replication 1-5 datasets
All six R scripts are available in the ./scripts/b_coh2_scripts/ inside the script dir, which can be run usinga  shell script
```
sbatch shell/coh2_block.sh
```
You can change the R script you call inside the shell script and run the complete model and replication dataset for the subset model
You can get the coheritability breakdown dataset for the four target traits in complete and replicated scenarios by running the .py script from the compute node in HPC 
```
python ./scripts/c_combine.py
```
## Pre-processing 
For this, you would need to get into the  compute node and an R session with R() in HPC and run the 
```
./scripts/c2_rerun_preprocessing.R
```
You will have all three wave-ratios selected and the necessary output for four target traits by running this prep-processing step to run the further step.
## Calculating BLUES and Heritability of synthetic trait/ selected wave ratios and each target trait in EF and MW locations
Again, you would need to run this in an interactive session in R
```
./scripts/d_rerun_modelfit.R
```
## Single-trait second stage model for all four target traits 
```
./scripts/e_run_singletraits_cv_completemodel.R
```
R script for calculating the accuracy and GebV of the single-trait for the complete model is listed above, using the below-listed bash script, which is modified for each target trait sla, narea, plsr-sla, plsr-narea, and can run the e_run_singletraits_cv_completemodel.R script
```
sbatch ./shell/st_cv_array.sh
```
## Multi-trait second stage model for all four target traits
R scripts for each trait can be found inside the scripts dir, and can be submitted through a bash script for each trait, which includes both CV1 and CV2 calculation
```
./scripts/f_completemodel_parallel/
```
bash script for submitting the R scripts above 
```
sbatch ./shell/run_stage2_narea.sh
```
## Selecting the synthetic trait with the lowest Coheritability 
Here again, for the pre-processing step and selecting the synthetic trait with the lowest coheritability, we need to run an interactive R session in HPC; it does need good memory and space because of all the clustering we do

```
./scripts/g_lowest_coh2/g_1_lowcoh2_preprocessingalltrait.R
```
Then, to fit the MT model with the lowest synthetic trait selected, you can submit the R script in the ./scripts/g_lowest_coh2/ dir using the below-mentioned bash scripts

```
sbatch ./shell/run_stage2_NT_lowcoh2.sh
```
There are other bash scripts already in the shell script that you can reuse if you dont like to change the wave_ratio value and R script to call everytime in bash script
```
sbatch ./shell/run_stage2_ST_lowcoh2.sh
sbatch ./shell/run_stage2_pnT_lowcoh2.sh
sbatch ./shell/run_stage2_psT_lowcoh2.sh
```
By this time, you should already have all the GEBV and accuracy files for the ST (single trait), MT(multi-trait), and S0 (MT model with synthetic trait's coh2=0)


## plotting 



