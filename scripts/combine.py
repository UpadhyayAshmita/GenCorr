# This file can be used to aggregate coheritability runs from multiple workers
import pandas as pd

trait = "ps"   # "narea", "sla", "pn", "ps"

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
