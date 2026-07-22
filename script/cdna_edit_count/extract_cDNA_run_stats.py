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
                                     [--dict_file <dict.csv.gz>]
                                     [--dict_type <regular|Dual_NEBuilder>]
                                     [--dict_name <RAD51C_XRCC3|...>]

If --out_csv is not supplied the CSV is written next to the job output file
with the same stem.

The Report_lines record which cDNA samples were parsed but not which dictionary
they were parsed against, so the dictionary is supplied here.  Pass --dict_file
and both the type and the name are read off its path; --dict_type/--dict_name
override that if the path is not in the canonical layout.
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

# Emitted only when step 01 ran with --cDNA_in_r2; the columns are added to the CSV
# only if at least one rep was actually screened, so unscreened batches are unchanged.
R2_COLUMNS = [
    "Read pairs screened on R2",
    "Read pairs rejected as dictionary contamination",
]

R2_FIELD_MAP = {
    "Reads_r2_screened":     "Read pairs screened on R2",
    "Reads_r2_dict_contam":  "Read pairs rejected as dictionary contamination",
}

FIELD_MAP = {
    "Reads_total":                              "Total reads",
    "Reads_barcode_matched":                    "Reads with barcode matched",
    "Reads_all_edit_pos_found":                 "Reads with template matched",
    "Reads_all_edit_pos_found_barcode_matched":  "Reads with barcode and template matched",
    "Barcodes":                                 "Unique barcodes in cDNA samples",
    "Barcodes_matched":                         "Unique barcodes matched in dictionary",
}


def dict_from_path(dict_file: str):
    """
    Split a dictionary path into (type, name).

    Canonical layout is
        .../barcode_to_variant/<libtype>/<form>/<name>.barcode_to_variant.csv.gz
    e.g. '.../barcode_to_variant/Dual_NEBuilder/concatenated/RAD51C_XRCC3.barcode_to_variant.csv.gz'
         -> ('Dual_NEBuilder', 'RAD51C_XRCC3')

    Returns (None, None) for the type / name that cannot be read off the path,
    leaving the caller free to fall back on an explicit flag.
    """
    name = os.path.basename(dict_file)
    for suffix in (".barcode_to_variant.csv.gz", ".barcode_to_variant.csv", ".csv.gz", ".csv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    # walk up past the form directory (concatenated/ combined/ blocks/) to the libtype
    parts = os.path.normpath(dict_file).split(os.sep)
    dtype = None
    if "barcode_to_variant" in parts:
        i = parts.index("barcode_to_variant")
        if i + 1 < len(parts) - 1:
            dtype = parts[i + 1]

    return dtype, (name or None)


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


def parse_job_output(job_output_file: str, dict_type=None, dict_name=None) -> pd.DataFrame:
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
            if dict_type is not None:
                row["Dictionary type"] = dict_type
            if dict_name is not None:
                row["Dictionary name"] = dict_name
            for field, col in FIELD_MAP.items():
                row[col] = raw.get(field)
            for field, col in R2_FIELD_MAP.items():
                row[col] = raw.get(field)

            records.append(row)

    if not records:
        sys.exit("No Report_line entries found in the job output file.")

    dict_cols = ([] if dict_type is None else ["Dictionary type"]) + \
                ([] if dict_name is None else ["Dictionary name"])
    # only surface the R2-screen columns when the screen actually ran
    screened = any(r.get("Read pairs screened on R2") for r in records)
    r2_cols = R2_COLUMNS if screened else []
    if not screened:
        for r in records:
            for col in R2_COLUMNS:
                r.pop(col, None)
    df = pd.DataFrame(records, columns=["Sample"] + dict_cols + COLUMNS + r2_cols)
    df = df.set_index("Sample")
    return df


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("job_output_file", help="Path to the job output file (e.g. muliPPIseq_DMS_editRate.o233697)")
    parser.add_argument("--out_csv",  default=None, help="Path for CSV output  (default: <job_output_stem>_stats.csv)")
    parser.add_argument("--dict_file", default=None,
                        help="Dictionary the samples were parsed against; type and name are read off its path")
    parser.add_argument("--dict_type", default=None,
                        help="Dictionary type (regular / Dual_NEBuilder); overrides --dict_file")
    parser.add_argument("--dict_name", default=None,
                        help="Dictionary name (RAD51C_XRCC3, RAD51D_XRCC2, ...); overrides --dict_file")
    args = parser.parse_args()

    stem = os.path.splitext(args.job_output_file)[0]
    out_csv  = args.out_csv  or f"{stem}_stats.csv"

    dict_type, dict_name = (None, None)
    if args.dict_file:
        dict_type, dict_name = dict_from_path(args.dict_file)
    dict_type = args.dict_type or dict_type
    dict_name = args.dict_name or dict_name

    if args.dict_file and (dict_type is None or dict_name is None):
        print(f"WARNING: could not read dictionary type/name off {args.dict_file}; "
              f"pass --dict_type/--dict_name explicitly", file=sys.stderr)

    df = parse_job_output(args.job_output_file, dict_type=dict_type, dict_name=dict_name)

    # ----- CSV -----
    df.to_csv(out_csv)
    print(f"CSV  written to: {out_csv}")

    # ----- preview -----
    print("\n--- Summary ---")
    print(df.to_string())


if __name__ == "__main__":
    main()
