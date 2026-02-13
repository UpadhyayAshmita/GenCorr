## How to reproduce the results
- High-performance computing with high memory and enough to submit the job as per the memory and resource allocated in a bash script in shell dir
- Also, some interactive R session is necessary to pre-process the data
## Clone the dir and create this dir structure 
- After cloning the dir, phenotypic data which is available in _ needs to be downloaded and put it inside data dir along with the genmoic data which is also available here _ needs to be downloaded.
- Dir structure should be like this:
```
GenCorr/
├── data/                    # phenotypic + genomic data
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



## For partitioning the complete vs replication 1-5 dataset and individual 
