#!/usr/bin/env python3
import os
import sys
import re
import gzip
import time
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length for Hamming distance")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Process PPIseq DMS dictionary build step 1 (v2: optional bc2 length modes)"
)
parser.add_argument(
    "--r1_fastq",
    required=True,
    help='Input R1 fastq file path. Input file prefix as in "<gene>_DMS_<block>_Dictionary". Example: RAD51D_DMS_1_Dictionary',
)
parser.add_argument(
    "--r2_fastq",
    required=True,
    help="Input R2 fastq file path",
)
parser.add_argument(
    "--r1_linker",
    required=True,
    help="Constant R1 sequence between bc1 and bc2 (split anchor)",
)
parser.add_argument(
    "--r2_fixed_seq_lib",
    required=True,
    help="File containing fixed sequence in R2",
)
parser.add_argument(
    "--r1_barcode2_lib",
    default=None,
    help="Optional file containing read 1 barcode 2 library (expected bc2 vs reference)",
)
parser.add_argument(
    "--bc1_min_len",
    type=int,
    default=12,
    help="Minimum length for barcode 1 (default: 12)",
)
parser.add_argument(
    "--bc1_max_len",
    type=int,
    default=15,
    help="Maximum length for barcode 1 (default: 15)",
)
parser.add_argument(
    "--bc2_set_len",
    type=int,
    default=None,
    help="Fixed length for barcode 2; not required if --r1_barcode2_lib matches gene/block (length taken from library)",
)
parser.add_argument(
    "--bc2_min_len",
    type=int,
    default=None,
    help="Minimum length for barcode 2; not required if --r1_barcode2_lib matches (see --bc2_max_len, --r1_post_bc2)",
)
parser.add_argument(
    "--bc2_max_len",
    type=int,
    default=None,
    help="Maximum length for barcode 2; not required if --r1_barcode2_lib matches",
)
parser.add_argument(
    "--r1_post_bc2",
    default=None,
    help="Constant R1 sequence immediately 3' of bc2; required with --bc2_min_len/--bc2_max_len to delimit bc2",
)
parser.add_argument(
    "--out",
    required=True,
    help="Output FASTA path; gzip-compressed if the path ends with .gz (case-insensitive)",
)

args = parser.parse_args()

has_range = args.bc2_min_len is not None and args.bc2_max_len is not None
has_partial_range = (args.bc2_min_len is not None) ^ (args.bc2_max_len is not None)
has_set = args.bc2_set_len is not None

if has_partial_range:
    logging.error(
        "Both --bc2_min_len and --bc2_max_len must be provided together, or use --bc2_set_len instead."
    )
    sys.exit(1)

# Construct input filenames from prefix
infile = args.r1_fastq
infile2 = args.r2_fastq
outfile = args.out
r2_fixed_seq_lib = args.r2_fixed_seq_lib
r1_barcode2_lib = args.r1_barcode2_lib

bc1_min_len = args.bc1_min_len
bc1_max_len = args.bc1_max_len
bc2_set_len = args.bc2_set_len
bc2_min_len = args.bc2_min_len
bc2_max_len = args.bc2_max_len

r1_linker = args.r1_linker.upper()
r1_post_bc2 = args.r1_post_bc2.upper() if args.r1_post_bc2 else None

# Extract gene name and block number from R1 filename.
# Parse the basename only: a directory in the path may contain underscores
# (e.g. .../Dual_NEBuilder/), which would otherwise truncate the gene token.
match = re.search(r"([^_]+)_DMS_(\d+)_Dictionary", os.path.basename(infile))
if match:
    gene = match.group(1)
    blk = match.group(2)
else:
    logging.error("Block number not found in R1 filename (expected pattern *_DMS_<n>_Dictionary).")
    sys.exit(1)

# Load r1 barcode2 sequence (needed before bc2 length mode can be resolved from library)
r1_bc2_seq = None
if r1_barcode2_lib is not None:
    with open(r1_barcode2_lib, "r") as f:
        next(f)  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            genename, block, fixed_seq = line.split(",")
            if (block == blk) and (genename == gene):
                r1_bc2_seq = fixed_seq.upper()

bc2_len_from_lib = len(r1_bc2_seq) if r1_bc2_seq is not None else None

if bc2_len_from_lib is not None:
    L = bc2_len_from_lib
    if args.bc2_set_len is not None and args.bc2_set_len != L:
        logging.warning(
            "--bc2_set_len (%s) does not match length of r1 barcode2 library sequence (%s); using library length.",
            args.bc2_set_len,
            L,
        )
    if has_range and not (bc2_min_len <= L <= bc2_max_len):
        logging.warning(
            "Barcode2 library sequence length (%s) is outside --bc2_min_len/--bc2_max_len (%s-%s); using library length for extraction.",
            L,
            bc2_min_len,
            bc2_max_len,
        )
    bc2_mode = "set"
    bc2_set_len = L
elif has_range and has_set:
    bc2_mode = "range"
elif has_range:
    bc2_mode = "range"
elif has_set:
    bc2_mode = "set"
else:
    logging.error(
        "Must provide either --r1_barcode2_lib with a matching gene/block row, "
        "or --bc2_set_len, or both --bc2_min_len and --bc2_max_len."
    )
    sys.exit(1)

if bc2_mode == "range":
    if args.r1_post_bc2 is None:
        logging.error(
            "--r1_post_bc2 is required when using --bc2_min_len and --bc2_max_len."
        )
        sys.exit(1)
    if args.bc2_min_len > args.bc2_max_len:
        logging.error(
            "--bc2_min_len (%s) cannot be greater than --bc2_max_len (%s).",
            args.bc2_min_len,
            args.bc2_max_len,
        )
        sys.exit(1)

# Load r2 fixed sequence
r2_fixed_seq = None
with open(r2_fixed_seq_lib, "r") as f:
    next(f)  # Skip header
    for line in f:
        line = line.strip()
        if not line:
            continue
        genename, block, fixed_seq = line.split(",")
        if (block == blk) and (genename == gene):
            r2_fixed_seq = fixed_seq.upper()

if r2_fixed_seq is None:
    logging.error("R2 fixed sequence not found for gene %s block %s", gene, blk)
    sys.exit(1)

logging.info("Gene: %s", gene)
logging.info("Block number: %s", blk)
logging.info("R1 linker (bc1|bc2): %s", r1_linker)
logging.info("BC2 mode: %s", bc2_mode)
if bc2_mode == "range":
    logging.info("R1 post-bc2: %s", r1_post_bc2)
    logging.info("BC2 length range: %s-%s", bc2_min_len, bc2_max_len)
else:
    logging.info("BC2 set length: %s", bc2_set_len)
logging.info("R2 fixed seq: %s", r2_fixed_seq)

# Open input files (handle .gz files)
if infile.endswith(".gz"):
    f_in = gzip.open(infile, "rt")
else:
    f_in = open(infile, "r")

if infile2.endswith(".gz"):
    f_in2 = gzip.open(infile2, "rt")
else:
    f_in2 = open(infile2, "r")

# Open output FASTA (gzip when filename ends with .gz)
if outfile.lower().endswith(".gz"):
    f_out = gzip.open(outfile, "wt")
else:
    f_out = open(outfile, "w")

total_ct = 0
pass_ct = 0
barcodes_found = 0
oligo_found = 0
start_time = time.time()

try:
    while True:
        # Read R1 (4 lines per read: ID, sequence, +, quality)
        id1 = f_in.readline()
        if not id1:
            break
        seq1 = f_in.readline()
        f_in.readline()  # Skip +
        f_in.readline()  # Skip quality

        # Read R2 (4 lines per read: ID, sequence, +, quality)
        id2 = f_in2.readline()
        if not id2:
            break
        seq2 = f_in2.readline()
        f_in2.readline()  # Skip +
        f_in2.readline()  # Skip quality

        total_ct += 1
        pass_flag = True
        bc1 = None
        bc2 = None
        oligo_seq = None

        seq1 = seq1.strip()
        # Split on linker (bc1 | linker | bc2 | …)
        parts = seq1.split(r1_linker)
        if len(parts) == 2:
            bc1 = parts[0]
            after_fixed = parts[1]

            if bc2_mode == "set":
                bc2 = after_fixed[:bc2_set_len]
            else:
                if r1_post_bc2 not in after_fixed:
                    pass_flag = False
                else:
                    head, _sep, _rest = after_fixed.partition(r1_post_bc2)
                    bc2 = head
                    bc2_len = len(bc2)
                    if not (bc2_min_len <= bc2_len <= bc2_max_len):
                        pass_flag = False

            if pass_flag and bc2 is not None:
                bc1_len = len(bc1)
                if (bc1_len >= bc1_min_len) and (bc1_len <= bc1_max_len):
                    if r1_bc2_seq is None:
                        barcodes_found += 1
                    elif len(bc2) == len(r1_bc2_seq) and hamming_distance(bc2, r1_bc2_seq) <= 1:
                        bc2 = r1_bc2_seq
                        barcodes_found += 1
                    else:
                        pass_flag = False
                else:
                    pass_flag = False
        else:
            pass_flag = False

        seq2 = seq2.strip()
        # Split on fixed sequence to extract oligo
        parts2 = seq2.split(r2_fixed_seq)
        if len(parts2) == 2:
            temp = parts2[0]
            oligo_seq = parts2[1]
            temp_len = len(temp)
            if temp_len > 0:
                oligo_found += 1
            else:
                pass_flag = False
        else:
            pass_flag = False

        if pass_flag and oligo_seq is not None:
            pass_ct += 1
            f_out.write(">" + bc1 + "_" + bc2 + "_" + str(pass_ct) + "\n")
            f_out.write(oligo_seq + "\n")

        if total_ct % 1000000 == 0:
            diff = time.time() - start_time
            print(
                "Processed reads: {}, passed reads: {}, time: {:.2f}".format(
                    total_ct, pass_ct, diff
                )
            )

finally:
    f_in.close()
    f_in2.close()
    f_out.close()

print("File_name\tTotal_reads\tbarcodes_found\toligo_found\tboth_found")
print("{}\t{}\t{}\t{}\t{}".format(infile, total_ct, barcodes_found, oligo_found, pass_ct))
