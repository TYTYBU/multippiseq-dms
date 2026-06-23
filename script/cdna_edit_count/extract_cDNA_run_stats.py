#!/usr/bin/env python3
"""
Parse cDNA run stats from a job output file and produce a summary CSV.

Each Report_line in the job output is tab-separated:
    Report_line  <filepath>  Reads_total  <N>  Reads_barcode_matched  <N>
                 Reads_all_edit_pos_found  <N>
                 Reads_all_edit_pos_found_barcode_matched  <N>
                 Barcodes  <N>  Barcodes_matched  <N>

Usage:
    python extract_cDNA_run_stats.py <job_output_file> [--out_csv <output.csv>]

If --out_csv is not supplied the CSV is written next to the job output file
with the same stem.
"""

import argparse
import os
import sys

import pandas as pd


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

COLUMNS = [
    "Total reads",
    "Reads with barcode matched",
    "Reads with template matched",
    "Reads with barcode and template matched",
    "Unique barcodes in cDNA samples",
    "Unique barcodes matched in dictionary",
]

FIELD_MAP = {
    "Reads_total":                              "Total reads",
    "Reads_barcode_matched":                    "Reads with barcode matched",
    "Reads_all_edit_pos_found":                 "Reads with template matched",
    "Reads_all_edit_pos_found_barcode_matched":  "Reads with barcode and template matched",
    "Barcodes":                                 "Unique barcodes in cDNA samples",
    "Barcodes_matched":                         "Unique barcodes matched in dictionary",
}


def sample_name_from_path(filepath: str) -> str:
    """
    Derive sample name from a cDNA_parsed output filepath.
    e.g. '.../NR-TV_HCT_RAD51C_XRCC3_Rep1_.cDNA_parsed.txt'
         -> 'NR-TV_HCT_RAD51C_XRCC3_Rep1'
    """
    basename = os.path.basename(filepath)
    name = basename.replace(".cDNA_parsed.csv", "").replace(".cDNA_parsed.txt", "")
    name = name.rstrip("_")
    return name


def parse_job_output(job_output_file: str) -> pd.DataFrame:
    records = []
    with open(job_output_file, "r") as fh:
        for line in fh:
            if not line.startswith("Report_line"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 14:
                print(f"WARNING: unexpected Report_line format, skipping: {line.strip()}", file=sys.stderr)
                continue

            filepath = parts[1]
            sample   = sample_name_from_path(filepath)

            # parts[2], parts[4], ... are field names; parts[3], parts[5], ... are values
            raw = {}
            for i in range(2, len(parts) - 1, 2):
                key = parts[i]
                try:
                    raw[key] = int(parts[i + 1])
                except (ValueError, IndexError):
                    raw[key] = None

            row = {"Sample": sample}
            for field, col in FIELD_MAP.items():
                row[col] = raw.get(field)

            records.append(row)

    if not records:
        sys.exit("No Report_line entries found in the job output file.")

    df = pd.DataFrame(records, columns=["Sample"] + COLUMNS)
    df = df.set_index("Sample")
    return df


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("job_output_file", help="Path to the job output file (e.g. muliPPIseq_DMS_editRate.o233697)")
    parser.add_argument("--out_csv",  default=None, help="Path for CSV output  (default: <job_output_stem>_stats.csv)")
    args = parser.parse_args()

    stem = os.path.splitext(args.job_output_file)[0]
    out_csv  = args.out_csv  or f"{stem}_stats.csv"

    df = parse_job_output(args.job_output_file)

    # ----- CSV -----
    df.to_csv(out_csv)
    print(f"CSV  written to: {out_csv}")

    # ----- preview -----
    print("\n--- Summary ---")
    print(df.to_string())


if __name__ == "__main__":
    main()
