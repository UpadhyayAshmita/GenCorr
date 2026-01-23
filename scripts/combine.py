#!/usr/bin/env python3
# Combine coheritability worker outputs without pandas/numpy (HPC-safe)

import os
import sys
import csv
import argparse
from glob import glob

COLUMNS = [
    "wave_1", "wave_2", "trait",
    "coh2", "h2_trait", "h2_ratio",
    "corg", "corgblup", "covs",
    "varw", "vars", "vartrait"
]

def die(msg: str, code: int = 1):
    print(msg, file=sys.stderr)
    sys.exit(code)

def iter_rows_whitespace(path: str):
    with open(path, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            yield line.split()

def iter_rows_csv(path: str):
    with open(path, "r", newline="") as fin:
        reader = csv.reader(fin)
        for row in reader:
            if not row:
                continue
            # handle accidental extra spaces around commas
            yield [x.strip() for x in row]

def combine_one_trait(trait: str, in_dir: str, out_file: str,
                      workers: int, chunk: int, wave_start: int,
                      columns, sep_mode: str, quiet: bool):
    if not os.path.isdir(in_dir):
        die(f"ERROR: input directory not found: {in_dir}")

    seen = set()
    kept = skipped_dups = total_rows = 0

    row_iter = iter_rows_whitespace if sep_mode == "whitespace" else iter_rows_csv

    with open(out_file, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(columns)

        for i in range(workers):
            start = (i * chunk) + wave_start
            end = ((i + 1) * chunk) + wave_start
            fn = os.path.join(in_dir, f"{trait}_{start}_{end}.csv")

            if not os.path.exists(fn):
                die(f"ERROR: missing worker file: {fn}")

            for parts in row_iter(fn):
                if len(parts) != len(columns):
                    die(
                        f"ERROR: unexpected column count in {fn}\n"
                        f"Got {len(parts)} fields, expected {len(columns)}\n"
                        f"Row: {parts}"
                    )

                total_rows += 1
                key = (parts[0], parts[1])  # wave_1, wave_2
                if key in seen:
                    skipped_dups += 1
                    continue
                seen.add(key)
                writer.writerow(parts)
                kept += 1

    if not quiet:
        print(f"Trait={trait}")
        print(f"Wrote: {out_file}")
        print(f"Total input rows: {total_rows}")
        print(f"Kept rows: {kept}")
        print(f"Skipped duplicates: {skipped_dups}")

        # print first 10 lines
        with open(out_file, "r") as f:
            for j, line in enumerate(f):
                print(line.rstrip("\n"))
                if j >= 10:
                    break

def main():
    ap = argparse.ArgumentParser(
        description="Combine GenCorr coheritability worker outputs (no pandas)."
    )
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--trait", help="Single trait to combine (e.g., narea)")
    g.add_argument("--traits", nargs="+", help="Multiple traits to combine (e.g., narea sla pn ps)")

    ap.add_argument("--in-dir", default="./output_complete", help="Input directory with worker files")
    ap.add_argument("--out-dir", default=".", help="Output directory")
    ap.add_argument("--workers", type=int, default=50, help="Number of worker files")
    ap.add_argument("--chunk", type=int, default=43, help="Wavelengths per worker")
    ap.add_argument("--wave-start", type=int, default=350, help="Starting wavelength index")
    ap.add_argument("--sep", choices=["whitespace", "comma"], default="whitespace",
                    help="Input file delimiter: whitespace (default) or comma")
    ap.add_argument("--quiet", action="store_true", help="Less printing")

    args = ap.parse_args()

    traits = [args.trait] if args.trait else args.traits
    os.makedirs(args.out_dir, exist_ok=True)

    for tr in traits:
        out_file = os.path.join(args.out_dir, f"{tr}_breakdown.csv")
        combine_one_trait(
            trait=tr,
            in_dir=args.in_dir,
            out_file=out_file,
            workers=args.workers,
            chunk=args.chunk,
            wave_start=args.wave_start,
            columns=COLUMNS,
            sep_mode=args.sep,
            quiet=args.quiet
        )

if __name__ == "__main__":
    main()
