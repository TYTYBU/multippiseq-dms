#!/usr/bin/env python3
"""DMS cDNA parse Step1 (v2): bc1/bc2 naming, csv dict, --ignore_bc1/--ignore_bc2. Legacy: multiPPIseq_DMS_cDNA_parse_Step1.py"""
import sys
import argparse
import logging
from itertools import product
import pandas as pd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Parse DMS cDNA data and generate alternative barcode files'
)
parser.add_argument(
    '--dict_in',
    required=True,
    help=(
        'Barcode dictionary: tab-separated .txt/.txt.gz, or comma-separated '
        '.csv/.csv.gz (e.g. barcode_to_variant from step3)'
    ),
)
parser.add_argument(
    '--edit_mode',
    required=True,
    choices=['A2G', 'T2C'],
    help='Edit mode: A2G or T2C'
)
parser.add_argument(
    '--out_prefix',
    required=True,
    help='Output file prefix for {out_prefix}.bc1_alts.txt and {out_prefix}.bc2_alts.txt'
)
ig = parser.add_mutually_exclusive_group()
ig.add_argument(
    '--ignore_bc1',
    action='store_true',
    help='Dictionary keys bc1 as ""; only expand bc2 barcodes (matches step3 --ignore_bc1)',
)
ig.add_argument(
    '--ignore_bc2',
    action='store_true',
    help='Dictionary keys bc2 as ""; only expand bc1 barcodes (matches step3 --ignore_bc2)',
)

args = parser.parse_args()

infile = args.dict_in
out_prefix = args.out_prefix
bc_edit_mode = args.edit_mode
ignore_bc1 = args.ignore_bc1
ignore_bc2 = args.ignore_bc2

# Set up edit mode parameters
if bc_edit_mode == "T2C":
    ref_n = "T"
    alt_n = "C"
elif bc_edit_mode == "A2G":
    ref_n = "A"
    alt_n = "G"
else:
    logger.error("Edit mode not recognized. Must be A2G or T2C")
    sys.exit(1)


def _dict_sep(path: str) -> str:
    """Comma for .csv / .csv.gz; tab for legacy .txt / .txt.gz / other."""
    pl = path.lower()
    if pl.endswith(".csv.gz") or pl.endswith(".csv"):
        return ","
    return "\t"


# Read input: compression='infer' handles .gz; sep follows extension
logger.info(f"Loading DMS dictionary {infile}...")
try:
    df_dict = pd.read_csv(
        infile,
        sep=_dict_sep(infile),
        header=None,
        compression="infer",
    )
except IOError as e:
    logger.error("can't open {}: {}".format(infile, e))
    sys.exit(1)
ncols = df_dict.shape[1]

# Extract bc1 and bc2 barcodes
if ncols == 5 and df_dict[0].astype(str).str.contains('_').all():
    df_dict.columns = ['barcode_pair', 'gene_aa_str', 'support_reads', 'total_reads', 'accuracy']
    pair_parts = df_dict['barcode_pair'].astype(str).str.split('_', n=1, expand=True)
    if ignore_bc1:
        bc1_set = {""}
        bc2_set = set(pair_parts[1].fillna("").astype(str))
    elif ignore_bc2:
        bc1_set = set(pair_parts[0].fillna("").astype(str))
        bc2_set = {""}
    else:
        bc1_set = set(pair_parts[0].fillna("").astype(str))
        bc2_set = set(pair_parts[1].fillna("").astype(str))
elif ncols == 7:
    df_dict.columns = ['bc1', 'bc1_orf', 'bc2', 'bc2_orf_aa_str', 'support_reads', 'total_reads', 'accuracy']
    if ignore_bc1:
        bc1_set = {""}
        bc2_set = set(df_dict['bc2'].fillna("").astype(str))
    elif ignore_bc2:
        bc1_set = set(df_dict['bc1'].fillna("").astype(str))
        bc2_set = {""}
    else:
        bc1_set = set(df_dict['bc1'].fillna("").astype(str))
        bc2_set = set(df_dict['bc2'].fillna("").astype(str))
else:
    raise ValueError(f"Unexpected file format: {ncols} columns detected. Expected 5 or 7 columns.")

if ignore_bc1:
    logger.info("bc1 ignored (--ignore_bc1): only bc2 alt expansion")
elif ignore_bc2:
    logger.info("bc2 ignored (--ignore_bc2): only bc1 alt expansion")

logger.info(
    f"Dictionary entries: {len(df_dict)}, unique bc1 barcodes: {len(bc1_set)}, "
    f"unique bc2 barcodes: {len(bc2_set)}"
)


def all_substitutions(seq, ref_n, alt_n, include_original=False):
    """
    For a given sequence, return the set of all unique sequences where any subset (including all)
    of the ref_n bases are changed to alt_n. If include_original is True, include the original sequence.
    If ref_n is not present, returns an empty set (or set containing the original if include_original).
    """
    # Find positions of ref_n in seq
    positions = [i for i, base in enumerate(seq) if base == ref_n]
    n = len(positions)
    out = set()

    if n == 0:
        if include_original:
            out.add(seq)
        return out

    # For each possible combination of positions
    # Each position can be swapped (True) or kept (False)
    for mask in product([False, True], repeat=n):
        if not any(mask) and not include_original:
            continue  # skip all ref_n kept (original sequence) unless include_original
        seq_list = list(seq)
        for j, swap in enumerate(mask):
            if swap:
                seq_list[positions[j]] = alt_n
        out.add(''.join(seq_list))
    return out


def _write_bc_alts(path: str, bc_list: list, ref_n: str, alt_n: str) -> None:
    """Append tab bc_alt, bc rows to path (same format as before)."""
    for i, bc in enumerate(bc_list):
        bc_alts = all_substitutions(bc, ref_n, alt_n, include_original=True)
        bc_alts = list(bc_alts)
        df_alts = pd.DataFrame({'bc_alt': bc_alts, 'bc': bc})
        mode = 'w' if i == 0 else 'a'
        df_alts.to_csv(path, mode=mode, header=False, index=False, sep='\t')


bc1_alts_path = out_prefix + ".bc1_alts.txt"
bc2_alts_path = out_prefix + ".bc2_alts.txt"

# Barcodes may contain edits as well — infer all possible alt barcodes from the dictionary
logger.info("Infering all possible bc1 barcodes...")
_write_bc_alts(bc1_alts_path, sorted(bc1_set), ref_n, alt_n)

logger.info("Infering all possible bc2 barcodes...")
_write_bc_alts(bc2_alts_path, sorted(bc2_set), ref_n, alt_n)
