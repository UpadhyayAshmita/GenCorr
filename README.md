## How to reproduce the results
- High-performance computing with high memory and enough to submit the job as per the memory and resource allocated in a bash script in the shell directory
- Also, some interactive R sessions are necessary to pre-process the data and get some output for serial computing in several steps
## Clone the dir and create this dir structure 
- After cloning the dir, phenotypic data, which is available in _ needs to be downloaded and put inside the data dir along with the genomic data, which is also available here _, and Names_WEST_SF.csv data for filtering later on the steps needs to be downloaded.
- Dir structure should be like this:
```
GenCorr/
├── data/                    # phenotypic + genomic data + Names_WEST.csv file 
├── figure/                  #for result figure 
├── function/                #custom functions
├── logs/                    #for checking job logs
├── output/                  #for output of pre-processing, first stage, and second stage model
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
source(./scripts/01_sampling_rep.R)
```
## Coheritability ratios generation for the complete dataset and replication 1-5 datasets
All six R scripts are available in the ./scripts/ inside the script dir, which can be run using shell script
```
source(./scripts/02_coh2_block.R)
source(./scripts/02_coh2_block1.R)
source(./scripts/02_coh2_block_rep2.R)
source(./scripts/02_coh2_block_rep3.R)
source(./scripts/02_coh2_block_rep4.R)
source(./scripts/02_coh2_block_rep5.R)
```
```
sbatch shell/coh2_block.sh  
```
You can change the R script you call inside the shell script and run the complete model and replication dataset for the subset model; just change the trait name and R script path in the bash script
You can get the coheritability breakdown dataset for the four target traits in complete and replicated scenarios by running the .py script from the compute node in HPC for complete model 
```
source(./scripts/03_combine.py)
# This file can be used to aggregate coheritability runs from multiple workers
import pandas as pd

trait = "sla"   # you can change the trait =  "narea", "sla", "pn", "ps" each time you run the script 

columns = ["wave_1", "wave_2", "trait",
           "coh2", "h2_trait", "h2_ratio",
           "corg", "corgblup", "covs",
           "varw", "vars", "vartrait"]

workers = 50

df_list = []
for i in range(workers):
    start = (i * 43) + 350
    end   = ((i + 1) * 43) + 350

    file_name = f"./output_complete/{trait}_{start}_{end}.csv"
    df = pd.read_csv(file_name, header=None, sep=" ", index_col=False, names=columns)
    df_list.append(df)

combined_df = pd.concat(df_list, ignore_index=True)
combined_df = combined_df.drop_duplicates(subset=["wave_1", "wave_2"], keep="first")
combined_df.to_csv(f"{trait}_breakdown.csv", index=False)

loaded_df = pd.read_csv(f"{trait}_breakdown.csv")
print(loaded_df.head(10))

```

But to get the breakdown data set for rep1 to rep 5, for eg narea 
```
# This file can be used to aggregate coheritability runs from multiple workers
import pandas as pd

trait = "narea"   # you can change the trait = to all four trait "narea", "sla", "pn", "ps"

columns = ["wave_1", "wave_2", "trait",
           "coh2", "h2_trait", "h2_ratio",
           "corg", "corgblup", "covs",
           "varw", "vars", "vartrait"]

workers = 50

df_list = []
for i in range(workers):
    start = (i * 43) + 350
    end   = ((i + 1) * 43) + 350

    file_name = f"./output1/{trait}_{start}_{end}.csv"   # remeber to direct the path to output dir where the coh2_block.R output were saved for rep1 for eg output1 for rep1 m output2 for rep2 likewise for original pipeline
    df = pd.read_csv(file_name, header=None, sep=" ", index_col=False, names=columns)
    df_list.append(df)

combined_df = pd.concat(df_list, ignore_index=True)
combined_df = combined_df.drop_duplicates(subset=["wave_1", "wave_2"], keep="first")
combined_df.to_csv(f"{trait}_breakdown1.csv", index=False) # remeber to name the file and  the path to trait_ breakdown1 for downstream pipeline
```
## Pre-processing 
For this, you would need to get into the  compute node and an R session with R() in HPC and run the script below:
```
source(./scripts/c2_rerun_preprocessing.R)
```
You will have all three wave-ratios selected and the necessary output for four target traits by running this prep-processing step to run the further step.
## Calculating BLUES and Heritability of synthetic trait/ selected wave ratios and each target trait in EF and MW locations
Again, you would need to run this in an interactive session in R
```
source(./scripts/d_rerun_modelfit.R)
```
## Single-trait second stage model for all four target traits 
```
source(./scripts/e_run_singletraits_cv_completemodel.R)
```
R script for calculating the accuracy and GebV of the single-trait for the complete model is listed above, using the below-listed bash script, which is modified for each target trait sla, narea, plsr-sla, plsr-narea, and can run the e_run_singletraits_cv_completemodel.R script
```
sbatch ./shell/st_cv_array.sh
```
## Multi-trait second stage model for all four target traits
R scripts for each trait can be found inside the scripts/f_completemodel_parallel/ dir, and can be submitted through a bash script for each trait, which includes both CV1 and CV2 calculation
```
source (./scripts/f_completemodel_parallel/rerun_narea_completemodel_efmw.R)
source (./scripts/f_completemodel_parallel/rerun_sla_completemodel_efmw.R)
source(./scripts/f_completemodel_parallel/rerun_ps_completemodel_efmw.R)
source( ./scripts/f_completemodel_parallel/rerun_ps_completemodel_efmw.R)
```
bash script for submitting the R scripts above 
```
sbatch ./shell/run_stage2_narea.sh
```
## Selecting the synthetic trait with the lowest Coheritability 
Here again, for the pre-processing step and selecting the synthetic trait with the lowest coheritability, we need to run an interactive R session in HPC; it does need good memory and space because of the clustering work we do to select the trait

```
./scripts/g_lowest_coh2/g_1_lowcoh2_preprocessingalltrait.R
```
Then, to fit the MT model with the lowest synthetic trait selected, you can submit the R script in the ./scripts/g_lowest_coh2/ dir using the below-mentioned bash scripts

```
sbatch ./shell/run_stage2_NT_lowcoh2.sh
```
There are other bash scripts already in the shell dir that you can reuse if you don't like to change the wave_ratio value and R script path to call every time in the  bash script
```
sbatch ./shell/run_stage2_ST_lowcoh2.sh
sbatch ./shell/run_stage2_pnT_lowcoh2.sh
sbatch ./shell/run_stage2_psT_lowcoh2.sh
```
By this time, you should already have all the GEBV and accuracy files for the ST (single trait), MT(multi-trait), and S0 (MT model with synthetic traits' coh2=0)
in the complete dataset model

## plotting 

Again, all the graphs and plots in the result section were generated using the R script below, which you can run interactively in R() on the compute node

```
./scripts/h_rerun_plot.R
```

### Replication 1-5 
Now, for selecting the three synthetic traits for each target trait in each rep, and to get BLUEs from the  first stage, and GEBV and accuracy from the second stage we have the following steps:

# Pre-processing

You need to call the function first to get the pre-processing output for all 1-5 reps, and you can run it in an interactive compute node in HPC
```
source(./function/aux_function.R)
nareapreprocessing_allrep (rep_id,  # while calling script change your rep_id=1 to 5 to get all rep pre-processing output
                    breakdown_dir = "./output",
                    phenotypes_path = "./data/phenotypes_whole.csv",
                    sample_path = "./data/sample_per_rep.csv",
                    out_dir = "./output",
                    k_groups = 3) #if you want to change your path and file name for input & output dir, change likewise
```

```
slapreprocessing_allrep(rep_id,
                        breakdown_dir = "./output",
                        phenotypes_path = "./data/phenotypes_whole.csv",
                        sample_path = "./data/sample_per_rep.csv",
                        out_dir = "./output",
                        k_groups = 3)
```
```
pnpreprocessing_allrep(rep_id,
                       breakdown_dir = "./output",
                       phenotypes_path = "./data/phenotypes_whole.csv",
                       sample_path = "./data/sample_per_rep.csv",
                       out_dir = "./output",
                       k_groups = 3,
                       pn_trait_col = "fs_plsr_narea")
```
```
pspreprocessing_allrep(rep_id,
                       breakdown_dir = "./output",
                       phenotypes_path = "./data/phenotypes_whole.csv",
                       sample_path = "./data/sample_per_rep.csv",
                       out_dir = "./output",
                       k_groups = 3,
                       ps_trait_col = "plsr_sla_sorghum")
```
## BLUEs /First stage model fitting for synthetic traits selected in four target traits over five reps
There are 3 synthetic traits selected for 4 targets in each rep, and we also have a synthetic trait with the lowest coheritability ~0 selected for each target trait in each rep. 
You can run the four scripts below to get the output/ Blues file for all the synthetic traits mentioned above :
```
source("./scripts/h_replication_alltraits_preprocesstostage1/narea_allrep_stage1.R)
source("./scripts/h_replication_alltraits_preprocesstostage1/sla_allrep_stage1.R)
source("./scripts/h_replication_alltraits_preprocesstostage1/ps_allrep_stage1.R)
source("./scripts/h_replication_alltraits_preprocesstostage1/pn_allrep_stage1.R)
```
## second stage model for replication 1-5 

To fit the single-trait (ST) model for all four target traits in 5 different reps, you can submit the R script for each trait; below are the R scripts:
```
source(./scripts/singletrait_narea_secondstageallrep.R)
source(./scripts/singletrait_sla_secondstageallrep.R)
source(./scripts/singletrait_pn_secondstageallrep.R)
source(./scripts/singletrait_ps_secondstageallrep.R)
```
You can call the R scripts above by calling the bash script from the shell dir 
```
source(./shell/singletrait_secondstageallrep.sh)
sbatch shell/singletrait_secondstageallrep.sh <targettrait> <rep> #target trait= narea, ps, pn, sla & rep= 1,2,3,4,5 
#sbatch shell/singletrait_secondstageallrep.sh narea 3  
#sbatch shell/singletrait_secondstageallrep.sh sla 5
#sbatch shell/singletrait_secondstageallrep.sh pn 2
#sbatch shell/singletrait_secondstageallrep.sh ps 1
```

Now to fit teh second stage model for all four target trait and their 3 synthetic trait each for 5 reps in two scenarios efmw and mwef, you can submit the R script below using the bash script:
```
source(./scripts/replication_secondstage/narea_combinedrep.R)
source(./scripts/replication_secondstage/sla_combinedrep.R)
source(./scripts/replication_secondstage/ps_combinedrep.R)
source(./scripts/replication_secondstage/pn_combinedrep.R)

```
You can submit the R script using a generic bash script for this step
```
source(./shell/secondstage_allrep.sh)
# Usage:
sbatch shell/secondstage_allrep.sh <target> <trait1> <trait2> <trait3> <rep> 
#eg: target= ps, pn, narea, sla & trait1=wave_1715_wave_1691, trait2= wave_728_wave_1071, trait3=wave_2335_wave_1446 & rep=3
```
