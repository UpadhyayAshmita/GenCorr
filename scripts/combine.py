#This file can be used to aggregate coheritability runs from multiple workers

import pandas as pd 
import math
import torch

# trait = "sla" # "narea", "sla", "pn", "ps"
trait = "narea"
# trait = "pn"
# trait = "ps"

traits = ["narea", "sla", "pn", "ps"]

columns = ["wave_1", "wave_2", "trait", 
"coh2", "h2_trait", "h2_ratio", 
"corg", "corgblup",  "covs",
"varw","vars", "vartrait"]
full_len = (2500-350) * 43

#number of files to read
workers = 50

for trait in traits:

    #empty tensor with -1 values
    data = torch.ones(((43 * workers)+1, 2500-350+1))# * -1  #plus one is to make sure 2500 and 350 are included!

    df_list = []
    #load data 
    for i in range(0,workers):
        start = (i * 43) + 350 
        end = ((i + 1) * 43) + 350
        shift = 350 #index shift

        file_name = f"./output_complete/{trait}_{start}_{end}.csv"
        df = pd.read_csv(file_name, header=None, sep=" ", index_col=False, names=columns)
        df_list.append(df)

    #combine dataframes
    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df = combined_df.drop_duplicates(subset=['wave_1', 'wave_2'], keep='first')

    ##check all values are there
    #for row in df.itertuples():
    #    wave_1 = int(row.wave_1[5:])
    #    wave_2 = int(row.wave_2[5:])
    #    data[wave_1-350, wave_2-350] = 1
    #print(sum(sum(data)))
    #print(data.shape[0] * data.shape[1])

    #save as csv file
    combined_df.to_csv(f'{trait}_breakdown.csv', index=False)

    #check loaded df
    loaded_df = pd.read_csv(f'{trait}_breakdown.csv')
    print(trait)
    print(loaded_df[:10])
