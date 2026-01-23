# scripts/combine.py
import pandas as pd
from pathlib import Path

# choose ONE trait to combine for now
trait = "narea"
traits = [trait]   # <-- key fix

columns = [
    "wave_1", "wave_2", "trait",
    "coh2", "h2_trait", "h2_ratio",
    "corg", "corgblup", "covs",
    "varw", "vars", "vartrait"
]

workers = 50
chunk = 43
wave_start = 350
in_dir = Path("output_complete")

for trait in traits:
    df_list = []

    for i in range(workers):
        start = (i * chunk) + wave_start
        end   = ((i + 1) * chunk) + wave_start

        file_path = in_dir / f"{trait}_{start}_{end}.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"Missing: {file_path}")

        # NOTE: keep sep=" " only if your files are space-delimited
        df = pd.read_csv(file_path, header=None, sep=" ", names=columns)
        df_list.append(df)

    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df = combined_df.drop_duplicates(subset=["wave_1", "wave_2"], keep="first")

    out_path = Path(f"{trait}_breakdown.csv")
    combined_df.to_csv(out_path, index=False)

    print(f"{trait}: wrote {out_path} with n={len(combined_df)} rows")
    print(combined_df.head(10))
