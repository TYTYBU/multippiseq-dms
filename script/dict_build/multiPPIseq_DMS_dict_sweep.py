#!/usr/bin/env python3
"""
Sweep step-05 thresholds (--min_support_reads x --min_accuracy) and record how
dictionary yield responds, to help choose parameters (esp. for Dual_NEBuilder).

Thresholds only affect the final filter, so for each block we load + redistribute
the barcode counts ONCE (reusing the step3 functions, so the result is identical
to the real pipeline), derive each barcode pair's post-redistribution
(best_count, accuracy) and pre-redistribution raw read total, then re-count the
whole (support x accuracy) grid in memory.

Per (block, min_support_reads, min_accuracy) row:
  variant_mapped_reads          reads assigned a variant (step-04 output; denom)
  reads_in_faithful_barcode_pairs  raw reads belonging to QC-passing pairs (numer)
  frac_reads_retained           reads_in_faithful_barcode_pairs / variant_mapped_reads
  observed_barcode_pairs        distinct bc1_bc2 pairs seen among those reads (denom)
  faithful_barcode_pairs           pairs passing support+accuracy (numer)
  frac_pairs_retained           faithful_barcode_pairs / observed_barcode_pairs

Notes:
- Valid = unique best variant (no tie) AND best_count >= S AND accuracy >= A,
  matching step3's filter_and_write.
- Support/accuracy use POST-redistribution counts (as the pipeline does);
  reads_in_faithful_barcode_pairs uses PRE-redistribution raw counts (real reads,
  not inflated by neighbor borrowing), consistent with dict_stats.
- --ignore_bc1/--ignore_bc2 must match how step 05 was run (default: --ignore_bc2)
  so the sweep reproduces the actual dictionary at the same thresholds.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import logging
import os
import sys
from typing import Dict, List, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

VM_SUFFIX = "_Dictionary.block_ref_seqs.csv.gz"

FIELDNAMES = [
    "block",
    "min_support_reads",
    "min_accuracy",
    "variant_mapped_reads",
    "reads_in_faithful_barcode_pairs",
    "frac_reads_retained",
    "observed_barcode_pairs",
    "faithful_barcode_pairs",
    "frac_pairs_retained",
]


def load_step3(script_dir: str):
    path = os.path.join(script_dir, "multiPPIseq_DMS_dict_build_step3.py")
    spec = importlib.util.spec_from_file_location("step3", path)
    mod = importlib.util.module_from_spec(spec)
    # Register before exec so step3's multiprocessing workers (pickled by
    # qualified name "step3.<fn>") can be resolved in worker processes.
    sys.modules["step3"] = mod
    spec.loader.exec_module(mod)
    return mod


def pair_stats(
    step3,
    variant_map_in: str,
    *,
    ignore_bc1: bool,
    ignore_bc2: bool,
    jobs: int,
) -> Tuple[List[Tuple[int, float, int]], int]:
    """
    Returns (pairs, variant_mapped_reads).

    pairs: one tuple per (bc1, bc2) barcode pair that has a unique best variant:
           (best_count_postredist, accuracy_postredist, raw_total_preredist).
    Tied pairs are dropped here (they can never pass the pipeline filter) but are
    still counted in observed_barcode_pairs via len-of-all-pairs below... see caller.
    """
    hash_bc, bc1_set, bc2_set, _total, valid_lines = step3.load_raw_dict(
        variant_map_in, ignore_bc1=ignore_bc1, ignore_bc2=ignore_bc2
    )
    nbr1 = {} if ignore_bc1 else step3.build_hd1_neighbors(bc1_set, jobs)
    nbr2 = {} if ignore_bc2 else step3.build_hd1_neighbors(bc2_set, jobs)
    accum = step3.redistribute_counts(hash_bc, nbr1, nbr2)

    pairs: List[Tuple[int, float, int]] = []
    observed = 0
    for bc1, bc2_dict in accum.items():
        for bc2, var_ct in bc2_dict.items():
            observed += 1
            best_ct = 0
            total_ct = 0
            for c in var_ct.values():
                total_ct += c
                if c > best_ct:
                    best_ct = c
            n_tie = sum(1 for c in var_ct.values() if c == best_ct)
            if n_tie != 1:
                continue  # never callable; excluded from valid, kept in observed
            acc = best_ct / total_ct if total_ct else 0.0
            raw_total = sum(hash_bc[bc1][bc2].values())  # pre-redistribution reads
            pairs.append((best_ct, acc, raw_total))
    return pairs, observed, valid_lines


def sweep_block(
    step3,
    block: str,
    variant_map_in: str,
    supports: List[int],
    accuracies: List[float],
    *,
    ignore_bc1: bool,
    ignore_bc2: bool,
    jobs: int,
) -> List[Dict[str, object]]:
    logging.info("Block %s: loading + redistributing once", block)
    pairs, observed, vm_reads = pair_stats(
        step3, variant_map_in, ignore_bc1=ignore_bc1, ignore_bc2=ignore_bc2, jobs=jobs
    )
    logging.info(
        "Block %s: variant_mapped_reads=%s observed_barcode_pairs=%s callable_pairs=%s",
        block,
        vm_reads,
        observed,
        len(pairs),
    )
    rows: List[Dict[str, object]] = []
    for S in supports:
        for A in accuracies:
            n_valid = 0
            reads_valid = 0
            for best_ct, acc, raw_total in pairs:
                if best_ct >= S and acc >= A:
                    n_valid += 1
                    reads_valid += raw_total
            rows.append(
                {
                    "block": block,
                    "min_support_reads": S,
                    "min_accuracy": A,
                    "variant_mapped_reads": vm_reads,
                    "reads_in_faithful_barcode_pairs": reads_valid,
                    "frac_reads_retained": (
                        reads_valid / vm_reads if vm_reads else 0.0
                    ),
                    "observed_barcode_pairs": observed,
                    "faithful_barcode_pairs": n_valid,
                    "frac_pairs_retained": (
                        n_valid / observed if observed else 0.0
                    ),
                }
            )
    return rows


def discover_blocks(variant_mapped_dir: str, only) -> List[Tuple[str, str]]:
    out = []
    for fn in sorted(os.listdir(variant_mapped_dir)):
        if not fn.endswith(VM_SUFFIX):
            continue
        block = fn[: -len(VM_SUFFIX)]
        if only and block not in only:
            continue
        out.append((block, os.path.join(variant_mapped_dir, fn)))
    return out


def parse_args():
    p = argparse.ArgumentParser(
        description="Sweep support/accuracy thresholds; record dictionary yield -> CSV."
    )
    p.add_argument(
        "--variant_mapped_dir",
        required=True,
        help="variant_mapped/<LIBTYPE>/blocks (step-04 per-read CSVs)",
    )
    p.add_argument("--out", required=True, help="Output sweep CSV path")
    p.add_argument(
        "--script_dir",
        default=os.path.dirname(os.path.abspath(__file__)),
        help="Dir holding multiPPIseq_DMS_dict_build_step3.py (default: this script's dir)",
    )
    p.add_argument(
        "--support_values",
        type=int,
        nargs="+",
        default=[1, 2, 3, 5, 10],
        help="min_support_reads values to sweep (default: 1 2 3 5 10)",
    )
    p.add_argument(
        "--accuracy_values",
        type=float,
        nargs="+",
        default=[0.6, 0.7, 0.8, 0.9, 1.0],
        help="min_accuracy values to sweep (default: 0.6 0.7 0.8 0.9 1.0)",
    )
    p.add_argument("--blocks", nargs="*", default=None, help="Subset of block keys")
    ig = p.add_mutually_exclusive_group()
    ig.add_argument("--ignore_bc1", action="store_true")
    ig.add_argument("--ignore_bc2", action="store_true")
    p.add_argument("--jobs", "-j", type=int, default=1, help="Workers for HD1 search")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not os.path.isdir(args.variant_mapped_dir):
        logging.error("Not a directory: %s", args.variant_mapped_dir)
        sys.exit(1)

    step3 = load_step3(args.script_dir)
    blocks = discover_blocks(args.variant_mapped_dir, args.blocks)
    if not blocks:
        logging.error("No variant_mapped blocks found.")
        sys.exit(1)
    logging.info(
        "Sweeping %s block(s) x %s support x %s accuracy (ignore_bc1=%s ignore_bc2=%s)",
        len(blocks),
        len(args.support_values),
        len(args.accuracy_values),
        args.ignore_bc1,
        args.ignore_bc2,
    )

    all_rows: List[Dict[str, object]] = []
    for block, vm_path in blocks:
        all_rows.extend(
            sweep_block(
                step3,
                block,
                vm_path,
                sorted(args.support_values),
                sorted(args.accuracy_values),
                ignore_bc1=args.ignore_bc1,
                ignore_bc2=args.ignore_bc2,
                jobs=max(1, args.jobs),
            )
        )

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        for r in all_rows:
            r = dict(r)
            r["frac_reads_retained"] = f"{r['frac_reads_retained']:.6g}"
            r["frac_pairs_retained"] = f"{r['frac_pairs_retained']:.6g}"
            writer.writerow(r)
    logging.info("Wrote %s rows -> %s", len(all_rows), args.out)


if __name__ == "__main__":
    main()
