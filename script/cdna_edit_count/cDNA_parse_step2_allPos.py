#!/usr/bin/env python3
"""DMS cDNA parse Step2 (v2): bc1/bc2, csv dict, --ignore_bc1/--ignore_bc2. Legacy: multiPPIseq_DMS_cDNA_parse_Step2_allPos.py"""
import sys
import argparse
import gzip
import logging
import re
import time
import pandas as pd
from itertools import islice

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Parse DMS cDNA data and analyze edits at all positions'
)
parser.add_argument(
    '--cDNA_in',
    required=True,
    help='cDNA R1.fastq file (can be gzipped)'
)
parser.add_argument(
    '--dict_in',
    required=True,
    help='Input dictionary file (can be gzipped)'
)
parser.add_argument(
    '--bc_alt_prefix',
    required=True,
    help='Prefix for *.bc1_alts.txt and *.bc2_alts.txt (from cDNA_parse_Step1_v2)',
)
parser.add_argument(
    '--out',
    required=True,
    help='Output file path'
)
ig = parser.add_mutually_exclusive_group()
ig.add_argument(
    '--ignore_bc1',
    action='store_true',
    help='Use bc1="" for lookups (must match Step1_v2 --ignore_bc1)',
)
ig.add_argument(
    '--ignore_bc2',
    action='store_true',
    help='Use bc2="" for lookups (must match Step1_v2 --ignore_bc2)',
)
parser.add_argument(
    '--edit_mode',
    required=True,
    choices=['A2G', 'T2C'],
    help='Edit mode: A2G or T2C'
)
parser.add_argument(
    '--template_seq',
    required=True,
    help='Template sequence'
)
parser.add_argument(
    '--regex_match_len',
    type=int,
    default=8,
    help='Length of template prefix/suffix used to build start/end regex anchors (default: 8)',
)
parser.add_argument(
    '--r1_post_bc2',
    default=None,
    help='Optional R1 sequence after bc2; first regex_match_len bases form an anchor regex. '
    'If set, bc2 is taken from the read between the end of cDNA and this anchor.',
)
parser.add_argument('--bc1_min_len', type=int, default=None, help='Minimum bc1 length (inclusive)')
parser.add_argument('--bc1_max_len', type=int, default=None, help='Maximum bc1 length (inclusive)')
parser.add_argument('--bc2_set_len', type=int, default=None, help='Exact bc2 length (mutually exclusive with --bc2_min_len/--bc2_max_len)')
parser.add_argument('--bc2_min_len', type=int, default=None, help='Minimum bc2 length; use with --bc2_max_len')
parser.add_argument('--bc2_max_len', type=int, default=None, help='Maximum bc2 length; use with --bc2_min_len')

args = parser.parse_args()
if (args.bc2_min_len is None) ^ (args.bc2_max_len is None):
    parser.error('--bc2_min_len and --bc2_max_len must be given together')
if args.bc2_set_len is not None and (args.bc2_min_len is not None or args.bc2_max_len is not None):
    parser.error('use either --bc2_set_len or --bc2_min_len/--bc2_max_len, not both')

infile2 = args.cDNA_in
infile = args.dict_in
bc_alt_prefix = args.bc_alt_prefix
out_file = args.out
bc_edit_mode = args.edit_mode
template = args.template_seq
regex_match_len = args.regex_match_len
r1_post_bc2 = args.r1_post_bc2
bc1_min_len = args.bc1_min_len
bc1_max_len = args.bc1_max_len
bc2_set_len = args.bc2_set_len
bc2_min_len = args.bc2_min_len
bc2_max_len = args.bc2_max_len
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

# convert template to regex pattern by replacing ref_n with dot
template = template.upper()
template_regex_str = template.replace(ref_n, ".")

# Start/end regex anchors from template; extracted cDNA must match both and have len == template_len
# Ensure head and tail substrings cannot overlap.
template_len = len(template_regex_str)
head_len = min(regex_match_len, template_len)
tail_len = min(regex_match_len, template_len - head_len)
if head_len <= 0 or tail_len <= 0:
    parser.error(
        f"regex_match_len too large for template: template_len={template_len}, head_len={head_len}, tail_len={tail_len}"
    )
head_regex_str = template_regex_str[:head_len]
tail_regex_str = template_regex_str[template_len - tail_len :]
head_regex = re.compile(head_regex_str, re.IGNORECASE)
tail_regex = re.compile(tail_regex_str, re.IGNORECASE)
if r1_post_bc2:
    r1_post_bc2 = r1_post_bc2.upper()
    r1_post_bc2_regex_str = r1_post_bc2.replace(ref_n, ".")
    post_anchor_len = min(regex_match_len, len(r1_post_bc2_regex_str))
    post_regex_str = (
        r1_post_bc2_regex_str[:post_anchor_len] 
        if len(r1_post_bc2_regex_str) >= post_anchor_len 
        else r1_post_bc2_regex_str
    )
    post_bc2_regex = re.compile(post_regex_str, re.IGNORECASE)
    logger.info(
        'r1_post_bc2 anchor: first %d bases (escaped regex len=%d)',
        post_anchor_len,
        len(post_regex_str),
    )
else:
    post_bc2_regex = None

logger.info(
    'template_len=%d, regex_match_len=%d, head_regex_len=%d, tail_regex_len=%d',
    template_len,
    regex_match_len,
    len(head_regex_str),
    len(tail_regex_str),
)

# Infer edit positions based on ref_n in template_seq
edit_pos = []
for i, base in enumerate(template):
    if base == ref_n:
        edit_pos.append(i + 1)  # 1-based position


def _dict_sep(path: str) -> str:
    pl = path.lower()
    if pl.endswith(".csv.gz") or pl.endswith(".csv"):
        return ","
    return "\t"


# Read dictionary file as pandas DataFrame
df_dict = pd.read_csv(
    infile,
    sep=_dict_sep(infile),
    header=None,
    compression="infer",
)
ncols = df_dict.shape[1]

# Extract bc1 and bc2 barcodes
if ncols == 5 and df_dict[0].astype(str).str.contains('_').all():
    df_dict.columns = ['barcode_pair', 'gene_aa_str', 'support_reads', 'total_reads', 'accuracy']
    pair_parts = df_dict['barcode_pair'].astype(str).str.split('_', n=1, expand=True)
    bc1_list = pair_parts[0].fillna("").astype(str).tolist()
    bc2_list = pair_parts[1].fillna("").astype(str).tolist()
    outcome_list = df_dict['gene_aa_str'].tolist()
elif ncols == 7:
    df_dict.columns = ['bc1', 'bc1_orf', 'bc2', 'bc2_orf_aa_str', 'support_reads', 'total_reads', 'accuracy']
    bc1_list = df_dict['bc1'].fillna("").astype(str).tolist()
    bc2_list = df_dict['bc2'].fillna("").astype(str).tolist()
    outcome_list = (df_dict['bc1_orf'] + "__" + df_dict['bc2_orf_aa_str']).tolist()
else:
    raise ValueError(f"Unexpected file format: {ncols} columns detected. Expected 5 or 7 columns.")

logger.info(
    f"Dictionary entries: {len(df_dict)}, bc1 barcodes: {len(bc1_list)}, bc2 barcodes: {len(bc2_list)}"
)

# initialize result_dict
result_dict = {outcome: {
    'total': 0,
    'any_edited_ct': 0,
    'is_edited': 0,
    'edit_pos': {pos: 0 for pos in edit_pos}
} for outcome in set(outcome_list)}

# derive bc2 barcode length (0 if all empty / ignore_bc2)
bc2_bc_len = len(bc2_list[0]) if bc2_list else 0
if not ignore_bc2:
    if not all(len(bc) == bc2_bc_len for bc in bc2_list):
        raise ValueError(f"bc2 barcodes have different lengths: {bc2_list}")
    logger.info(f"bc2 barcode length: {bc2_bc_len}")
else:
    logger.info("bc2 ignored: bc2 barcode length treated as 0")


def bc1_len_ok(L):
    if ignore_bc1:
        return True
    if bc1_min_len is not None and L < bc1_min_len:
        return False
    if bc1_max_len is not None and L > bc1_max_len:
        return False
    return True


def bc2_len_ok(L):
    if ignore_bc2:
        return True
    if bc2_set_len is not None:
        return L == bc2_set_len
    if bc2_min_len is not None:
        return bc2_min_len <= L <= bc2_max_len
    # Default behavior: enforce dictionary-derived bc2 length.
    return L == bc2_bc_len


# create barcode pair dictionary from bc1 and bc2 barcodes and outcomes.
# A (bc1, bc2) pair can map to MULTIPLE outcomes when the dictionary pools
# per-block dicts (a bc1 reused across blocks for different variants), so map
# each pair to the SET of all its outcomes — a matching read is counted toward
# every variant involved.
barcode_pair_dict = {}
for (b1, b2, outcome) in zip(bc1_list, bc2_list, outcome_list):
    barcode_pair_dict.setdefault((b1, b2), set()).add(outcome)

# Pre-compute reverse lookup: for each bc1, which bc2 values are valid
bc1_to_bc2s = {}
for b1, b2 in zip(bc1_list, bc2_list):
    if b1 not in bc1_to_bc2s:
        bc1_to_bc2s[b1] = set()
    bc1_to_bc2s[b1].add(b2)

bc1_alts_path = bc_alt_prefix + ".bc1_alts.txt"
bc2_alts_path = bc_alt_prefix + ".bc2_alts.txt"

# Load all possible bc1 alt barcodes
try:
    f_in = open(bc1_alts_path, 'r')
except IOError as e:
    logger.error("can't open {}: {}".format(bc1_alts_path, e))
    sys.exit(1)

logger.info("Loading all possible bc1 alt barcodes...")
start_time = time.time()
total_ct = 0

bc1_alts_dict = {}
for line in f_in:
    line = line.rstrip('\n')
    if not line:
        continue
    parts = line.split('\t')
    if len(parts) >= 2:
        bc_alt = parts[0]
        bc = parts[1]
        if bc_alt not in bc1_alts_dict:
            bc1_alts_dict[bc_alt] = set()
        bc1_alts_dict[bc_alt].add(bc)
    
    total_ct += 1
    if total_ct % 100000 == 0:
        diff = time.time() - start_time
        logger.info("bc1 barcodes processed: {}, time: {:.2f}".format(total_ct, diff))
        start_time = time.time()

f_in.close()

# Load all possible bc2 alt barcodes
try:
    f_in = open(bc2_alts_path, 'r')
except IOError as e:
    logger.error("can't open {}: {}".format(bc2_alts_path, e))
    sys.exit(1)

logger.info("Loading all possible bc2 alt barcodes...")
start_time = time.time()
total_ct = 0

bc2_alts_dict = {}
for line in f_in:
    line = line.rstrip('\n')
    if not line:
        continue
    parts = line.split('\t')
    if len(parts) >= 2:
        bc_alt = parts[0]
        bc = parts[1]
        if bc_alt not in bc2_alts_dict:
            bc2_alts_dict[bc_alt] = set()
        bc2_alts_dict[bc_alt].add(bc)
    
    total_ct += 1
    if total_ct % 100000 == 0:
        diff = time.time() - start_time
        logger.info("bc2 barcodes processed: {}, time: {:.2f}".format(total_ct, diff))
        start_time = time.time()

f_in.close()

total_bc1_alts = len(bc1_alts_dict)
logger.info("Total bc1 alt barcodes: {}".format(total_bc1_alts))
total_bc2_alts = len(bc2_alts_dict)
logger.info("Total bc2 alt barcodes: {}".format(total_bc2_alts))

# Parse cDNA FASTQ file
if infile2.endswith('.gz'):
    try:
        f_in2 = gzip.open(infile2, 'rt')
    except IOError as e:
        logger.error("can't open pipe to {}: {}".format(infile2, e))
        sys.exit(1)
else:
    try:
        f_in2 = open(infile2, 'r')
    except IOError as e:
        logger.error("can't open {}: {}".format(infile2, e))
        sys.exit(1)

total_reads = 0
reads_cDNA_pass = 0
reads_barcode_found = 0
reads_cDNA_pass_barcode_found = 0
read_dict = {}  # read_dict[bc1_bc2] = barcode_pass (0 or 1)

logger.info("Parsing cDNA data...")
start_time = time.time()

# Use iterator for faster FASTQ reading
while True:
    # Read 4 lines at a time (FASTQ format)
    lines = list(islice(f_in2, 4))
    if not lines:
        break
    if len(lines) < 4:
        break
    seq = lines[1].rstrip('\n')
    total_reads += 1
    
    # Anchor cDNA with template start/end regex; slice length equals template_len
    match = head_regex.search(seq)
    if match:
        s = match.start()
        if s + template_len > len(seq):
            continue
        cDNA = seq[s : s + template_len]
        cDNA_len = len(cDNA)
        if cDNA_len != template_len:
            continue
        tail_seg = cDNA[-len(tail_regex_str) :]
        if not tail_regex.match(tail_seg):
            continue

        bc1_read = "" if ignore_bc1 else seq[:s]
        if not bc1_len_ok(len(bc1_read)):
            continue

        if ignore_bc2:
            bc2_read = ""
        elif r1_post_bc2:
            m_post = post_bc2_regex.search(seq, s + template_len)
            if not m_post:
                continue
            bc2_read = seq[s + template_len : m_post.start()]
        else:
            bc2_start = s + template_len
            bc2_end = bc2_start + bc2_bc_len
            if bc2_end > len(seq):
                continue
            bc2_read = seq[bc2_start:bc2_end]

        if not bc2_len_ok(len(bc2_read)):
            continue
        
        edit_ct = 0
        edit_dict = {pos: 0 for pos in edit_pos}
        cDNA_pass = 1
        
        for pos in edit_pos:
            if pos > cDNA_len:
                break  # do not consider positions over cDNA length
            
            n = cDNA[pos - 1]  # 1-based to 0-based conversion
            
            if (n != ref_n) and (n != alt_n):
                cDNA_pass = 0
                break
            
            if n == alt_n:
                edit_ct += 1
                edit_dict[pos] = 1
        
        if cDNA_pass:
            reads_cDNA_pass += 1
        
        barcode_pass = 0
        outcome_set = set()
        if (bc1_read in bc1_alts_dict) and (bc2_read in bc2_alts_dict):
            bc1_originals = bc1_alts_dict[bc1_read]
            bc2_originals = bc2_alts_dict[bc2_read]
            for original_bc1 in bc1_originals:
                if original_bc1 in bc1_to_bc2s:
                    valid_bc2s = bc1_to_bc2s[original_bc1]
                    matching_bc2s = bc2_originals & valid_bc2s
                    for original_bc2 in matching_bc2s:
                        outcome_set.update(barcode_pair_dict[(original_bc1, original_bc2)])
                        barcode_pass = 1
                        reads_barcode_found += 1
        
        bc_key = f"{bc1_read}_{bc2_read}"
        # Use max() to preserve barcode_pass=1 if any read with this barcode pair matched
        read_dict[bc_key] = max(read_dict.get(bc_key, 0), barcode_pass)
        
        if (cDNA_pass == 1) and (barcode_pass == 1):
            reads_cDNA_pass_barcode_found += 1
            
            for outcome in outcome_set:
                result_dict[outcome]['total'] += 1
                result_dict[outcome]['any_edited_ct'] += edit_ct
                if edit_ct > 0:
                    result_dict[outcome]['is_edited'] += 1
                
                for pos in edit_dict:
                    result_dict[outcome]['edit_pos'][pos] += edit_dict[pos]
    
    if total_reads % 1000000 == 0:
        diff = time.time() - start_time
        logger.info("cDNA reads parsed: {}, time: {:.2f}".format(total_reads, diff))
        start_time = time.time()

f_in2.close()

# Write output file
try:
    f_out = open(out_file, 'w')
except IOError as e:
    logger.error("can't open output file for writing: {}".format(e))
    sys.exit(1)

f_out.write("Outcome,Total_reads,Is_edited,Any_edit_ct")
for pos in edit_pos:
    f_out.write(f",Pos{pos}")
f_out.write("\n")

for outcome in sorted(result_dict.keys()):
    f_out.write(f"{outcome},{result_dict[outcome]['total']},{result_dict[outcome]['is_edited']},{result_dict[outcome]['any_edited_ct']}")

    for pos in edit_pos:
        f_out.write(f",{result_dict[outcome]['edit_pos'].get(pos, 0)}")
    f_out.write("\n")

f_out.close()

# Calculate statistics
read_bc_pair = len(read_dict)
read_bc_pair_in_dict = sum(read_dict.values())

# cDNA analysis report line
print("Report_line\t{}\tReads_total\t{}\tReads_barcode_matched\t{}\tReads_all_edit_pos_found\t{}\tReads_all_edit_pos_found_barcode_matched\t{}\tBarcodes\t{}\tBarcodes_matched\t{}".format(
    out_file,
    total_reads,
    reads_barcode_found,
    reads_cDNA_pass,
    reads_cDNA_pass_barcode_found,
    read_bc_pair,
    read_bc_pair_in_dict
))
