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
R script for calculating the accuracy and GebV of the single-trait for the complete model is listed above, using the below-listed bash script, which is modified for each target trait sla, narea, plsr-sla, plsr-narea can run the e_run_singletraits_cv_completemodel.R script
```
sbatch ./shell/st_cv_array.sh
```

