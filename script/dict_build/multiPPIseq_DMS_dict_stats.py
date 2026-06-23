#!/usr/bin/env python3
"""
Dictionary-build diagnostic stats (Python port of DMS_dict_stats.pl).

For each per-block dictionary, count how reads survive each stage of the
pipeline and how many high-confidence barcode pairs were called:

  total_reads                    raw read pairs in the R1 fastq
  reads_with_barcode             reads whose barcodes were extracted (step 02,
                                 processed_fa FASTA records)
  reads_with_barcode_and_variant reads also mapped to a known variant (step 04,
                                 variant_mapped CSV rows)
  reads_in_faithful_barcode_pairs   raw reads belonging to the barcode pairs that
                                 passed the step-05 QC filters (support/accuracy)
                                 -- i.e. the reads "captured" by the final dict
  faithful_barcode_pairs            barcode pairs in the final dictionary

NOTE on the renamed column: the Perl called `reads_in_faithful_barcode_pairs`
"treads_passed_QC", which was ambiguous -- it is not a per-read QC pass count
but the sum of raw read support over the barcode pairs that made it into the
final dictionary. (The Perl also *overwrote* rather than summed a pair's count
when the pair spanned multiple variants in the raw dict; this port sums them.)

Outputs a single CSV with one row per block. Blocks are discovered from the raw
R1 fastqs; the matching files in the other stage directories are located by the
`<gene>_DMS_<block>` key.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import multiprocessing as mp
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

FIELDNAMES = [
    "block",
    "total_reads",
    "reads_with_barcode",
    "reads_with_barcode_and_variant",
    "reads_in_faithful_barcode_pairs",
    "faithful_barcode_pairs",
]

R1_SUFFIX = "_Dictionary_R1.fastq.gz"


def _open_binary(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def count_lines(path: str) -> int:
    """Fast newline count (handles gzip). Counts a trailing unterminated line."""
    n = 0
    last = b"\n"
    with _open_binary(path) as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            n += chunk.count(b"\n")
            last = chunk[-1:]
    if last not in (b"\n", b""):
        n += 1
    return n


def open_text(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def reads_in_valid_pairs(raw_dict_path: str, final_dict_path: str) -> Tuple[int, int]:
    """
    Sum raw read support over the barcode pairs present in the final dictionary.

    raw_dict CSV  : bc1, bc2, variant, read_count   (one row per bc-pair x variant)
    final dict CSV: bc_pair (bc1_bc2), best_variant, best_ct, total_ct, accuracy

    Returns (reads_in_faithful_barcode_pairs, faithful_barcode_pairs).
    """
    pair_total: Dict[str, int] = {}
    with open_text(raw_dict_path) as fh:
        for row in csv.reader(fh):
            if len(row) < 4:
                continue
            bc1, bc2, _variant, ct = row[0], row[1], row[2], row[3]
            try:
                c = int(ct)
            except ValueError:
                continue
            pair_total[f"{bc1}_{bc2}"] = pair_total.get(f"{bc1}_{bc2}", 0) + c

    reads = 0
    valid_pairs = 0
    with open_text(final_dict_path) as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            bc_pair = row[0]
            valid_pairs += 1
            reads += pair_total.get(bc_pair, 0)
    return reads, valid_pairs


def stats_for_block(task: dict) -> Dict[str, object]:
    block = task["block"]
    logging.info("Processing block %s", block)
    total_reads = count_lines(task["raw_r1"]) // 4
    reads_with_bc = count_lines(task["processed_fa"]) // 2  # 2-line FASTA records
    reads_with_bc_var = count_lines(task["variant_mapped"])
    reads_valid, valid_pairs = reads_in_valid_pairs(
        task["raw_dict"], task["final_dict"]
    )
    return {
        "block": block,
        "total_reads": total_reads,
        "reads_with_barcode": reads_with_bc,
        "reads_with_barcode_and_variant": reads_with_bc_var,
        "reads_in_faithful_barcode_pairs": reads_valid,
        "faithful_barcode_pairs": valid_pairs,
    }


def discover_blocks(
    raw_dir: str,
    processed_fa_dir: str,
    variant_mapped_dir: str,
    blocks_dir: str,
    only: Optional[List[str]],
) -> List[dict]:
    tasks: List[dict] = []
    r1_files = sorted(
        f for f in os.listdir(raw_dir) if f.endswith(R1_SUFFIX)
    )
    for fn in r1_files:
        block = fn[: -len(R1_SUFFIX)]  # <gene>_DMS_<block>
        if only and block not in only:
            continue
        paths = {
            "block": block,
            "raw_r1": os.path.join(raw_dir, fn),
            "processed_fa": os.path.join(
                processed_fa_dir, f"{block}_Dictionary.fasta.gz"
            ),
            "variant_mapped": os.path.join(
                variant_mapped_dir, f"{block}_Dictionary.block_ref_seqs.csv.gz"
            ),
            "raw_dict": os.path.join(blocks_dir, f"{block}.raw_dict.csv.gz"),
            "final_dict": os.path.join(
                blocks_dir, f"{block}.barcode_to_variant.csv.gz"
            ),
        }
        missing = [
            k
            for k in ("processed_fa", "variant_mapped", "raw_dict", "final_dict")
            if not os.path.isfile(paths[k])
        ]
        if missing:
            logging.warning(
                "Skipping block %s; missing inputs: %s",
                block,
                ", ".join(paths[k] for k in missing),
            )
            continue
        tasks.append(paths)
    return tasks


def natural_key(block: str):
    m = re.match(r"(.+)_DMS_(\d+)$", block)
    if m:
        return (m.group(1), int(m.group(2)))
    return (block, 0)


def parse_args():
    p = argparse.ArgumentParser(
        description="Per-block dictionary-build diagnostic stats -> single CSV."
    )
    p.add_argument("--raw_dir", required=True, help="raw_data/<LIBTYPE> (R1 fastq.gz)")
    p.add_argument(
        "--processed_fa_dir", required=True, help="processed_fa/<LIBTYPE>"
    )
    p.add_argument(
        "--variant_mapped_dir",
        required=True,
        help="variant_mapped/<LIBTYPE>/blocks",
    )
    p.add_argument(
        "--blocks_dir",
        required=True,
        help="barcode_to_variant/<LIBTYPE>/blocks (holds *.raw_dict.csv.gz and *.barcode_to_variant.csv.gz)",
    )
    p.add_argument("--out", required=True, help="Output CSV path")
    p.add_argument(
        "--blocks",
        nargs="*",
        default=None,
        help="Optional subset of block keys (e.g. RAD51C_DMS_1); default: all found",
    )
    p.add_argument(
        "--jobs", "-j", type=int, default=1, help="Parallel workers over blocks"
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    for d in (
        args.raw_dir,
        args.processed_fa_dir,
        args.variant_mapped_dir,
        args.blocks_dir,
    ):
        if not os.path.isdir(d):
            logging.error("Not a directory: %s", d)
            sys.exit(1)

    tasks = discover_blocks(
        args.raw_dir,
        args.processed_fa_dir,
        args.variant_mapped_dir,
        args.blocks_dir,
        args.blocks,
    )
    if not tasks:
        logging.error("No blocks with complete inputs found.")
        sys.exit(1)
    logging.info("Computing stats for %s block(s)", len(tasks))

    jobs = max(1, args.jobs)
    if jobs == 1:
        rows = [stats_for_block(t) for t in tasks]
    else:
        with mp.Pool(min(jobs, len(tasks))) as pool:
            rows = pool.map(stats_for_block, tasks)

    rows.sort(key=lambda r: natural_key(r["block"]))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows -> %s", len(rows), args.out)


if __name__ == "__main__":
    main()
