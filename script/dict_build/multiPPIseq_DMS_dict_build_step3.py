#!/usr/bin/env python3
"""
Merged pipeline: former step3a + step3b + step3c (Perl).

Reads variant-mapped CSV (from step2), builds raw barcode×variant counts,
redistributes counts across Hamming-distance-1 barcode neighbors for bc1 and bc2
(unless --ignore_bc1 or --ignore_bc2; mutually exclusive). When ignored, that
coordinate is keyed as "". Then outputs high-fidelity barcode pairs passing
support and accuracy thresholds.

Read id in column 1 is parsed as bc1_bc2 or bc1_bc2_<suffix> (underscore-separated);
e.g. AGTGGGGCATCA_TTCAA_42 → bc1=AGTGGGGCATCA, bc2=TTCAA.

The original Perl step3c loaded the second barcode's neighbors into the wrong hash;
this script keeps bc1 and bc2 neighbor maps separate.
"""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import logging
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict
from typing import DefaultDict, Dict, Iterable, List, Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

BASES = ("A", "C", "G", "T")


def open_text(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def open_out(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def write_raw_dict_csv(
    hash_bc: Dict[str, Dict[str, Dict[str, int]]],
    path: str,
) -> None:
    """Pre-redistribution counts as CSV: bc1, bc2, variant, read_count (gzip if path ends with .gz, e.g. .csv.gz)."""
    logging.info("Writing raw dictionary CSV: %s", path)
    with open_out(path) as fh:
        writer = csv.writer(fh, delimiter=",", lineterminator="\n")
        for bc1 in sorted(hash_bc.keys()):
            for bc2 in sorted(hash_bc[bc1].keys()):
                for variant in sorted(hash_bc[bc1][bc2].keys()):
                    ct = hash_bc[bc1][bc2][variant]
                    writer.writerow([bc1, bc2, variant, ct])


def parse_args():
    p = argparse.ArgumentParser(
        description="Build barcode→variant dictionary with HD1 neighbor redistribution (step 3 merged)."
    )
    p.add_argument(
        "--variant_map_in",
        required=True,
        metavar="PATH",
        help="Variant map CSV: read_id (bc1_bc2[_suffix]), variant; gzip ok",
    )
    p.add_argument(
        "--out",
        "-o",
        required=True,
        help="Output filtered barcode→variant table as CSV (gzip if path ends with .gz)",
    )
    p.add_argument(
        "--min_support_reads",
        type=int,
        default=3,
        help="Minimum read support for best variant (default: 3)",
    )
    p.add_argument(
        "--min_accuracy",
        type=float,
        default=0.9,
        help="Minimum best_variant_count / total_count (default: 0.9)",
    )
    p.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=1,
        help="CPU workers for HD1 neighbor search (default: 1, single core)",
    )
    ig = p.add_mutually_exclusive_group()
    ig.add_argument(
        "--ignore_bc1",
        action="store_true",
        help="Key all rows with bc1=\"\"; skip bc1 HD1 neighbors and redistribution along bc1",
    )
    ig.add_argument(
        "--ignore_bc2",
        action="store_true",
        help="Key all rows with bc2=\"\"; skip bc2 HD1 neighbors and redistribution along bc2",
    )
    p.add_argument(
        "--raw_dict_out",
        metavar="PATH",
        default=None,
        help="Optional raw dict CSV: bc1, bc2, variant, read_count. Use .csv.gz (or any path ending in .gz) for gzip",
    )
    return p.parse_args()


def load_raw_dict(
    path: str,
    *,
    ignore_bc1: bool = False,
    ignore_bc2: bool = False,
) -> Tuple[
    Dict[str, Dict[str, Dict[str, int]]],
    Set[str],
    Set[str],
    int,
    int,
]:
    """
    Step 3a equivalent: aggregate (bc1, bc2) -> variant -> count.
    Each CSV line: first field = read id with bc1_bc2 or bc1_bc2_<suffix>, second = variant.
    If ignore_bc1, bc1 is stored as "". If ignore_bc2, bc2 is stored as "".
    """
    hash_bc: DefaultDict[str, DefaultDict[str, DefaultDict[str, int]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(int))
    )
    bc1_set: Set[str] = set()
    bc2_set: Set[str] = set()
    total_lines = 0
    valid_lines = 0
    t0 = time.time()
    with open_text(path) as f:
        reader = csv.reader(f)
        for row in reader:
            total_lines += 1
            if len(row) < 2:
                continue
            temp, variant = row[0], row[1]
            parts = temp.split("_")
            if len(parts) < 2:
                continue
            bc1, bc2 = parts[0], parts[1]
            if ignore_bc1:
                bc1 = ""
            if ignore_bc2:
                bc2 = ""
            bc1_set.add(bc1)
            bc2_set.add(bc2)
            hash_bc[bc1][bc2][variant] += 1
            valid_lines += 1
            if total_lines % 1_000_000 == 0:
                logging.info(
                    "Loaded CSV lines: %s (valid %s), %.2fs",
                    total_lines,
                    valid_lines,
                    time.time() - t0,
                )
                t0 = time.time()
    logging.info(
        "CSV total lines: %s, valid: %s, unique bc1: %s, unique bc2: %s",
        total_lines,
        valid_lines,
        len(bc1_set),
        len(bc2_set),
    )
    return hash_bc, bc1_set, bc2_set, total_lines, valid_lines


def _neighbors_one_length(args: Tuple[int, Tuple[str, ...]]) -> Dict[str, Set[str]]:
    """Worker: HD1 neighbors for every barcode in a length bucket."""
    _length, barcodes = args
    bset = set(barcodes)
    out: Dict[str, Set[str]] = {}
    for bc in barcodes:
        nbrs: Set[str] = set()
        for i, ch in enumerate(bc):
            for alt in BASES:
                if alt == ch:
                    continue
                nb = bc[:i] + alt + bc[i + 1 :]
                if nb in bset:
                    nbrs.add(nb)
        out[bc] = nbrs
    return out


def _neighbors_chunk_subset(args: Tuple[frozenset, Tuple[str, ...]]) -> Dict[str, Set[str]]:
    """Worker: HD1 neighbors for a subset of barcodes; membership uses full set for that length."""
    full_set, chunk = args
    out: Dict[str, Set[str]] = {}
    for bc in chunk:
        nbrs: Set[str] = set()
        for i, ch in enumerate(bc):
            for alt in BASES:
                if alt == ch:
                    continue
                nb = bc[:i] + alt + bc[i + 1 :]
                if nb in full_set:
                    nbrs.add(nb)
        out[bc] = nbrs
    return out


def build_hd1_neighbors(
    barcodes: Iterable[str], jobs: int
) -> Dict[str, Set[str]]:
    """Step 3b equivalent without partitions: map each barcode -> set of HD1 neighbors in the full set."""
    by_len: DefaultDict[int, List[str]] = defaultdict(list)
    for bc in barcodes:
        by_len[len(bc)].append(bc)

    merged: Dict[str, Set[str]] = {}
    if not by_len:
        return merged

    jobs = max(1, jobs)
    chunk_tasks: List[Tuple[frozenset, Tuple[str, ...]]] = []

    for L, seqs in by_len.items():
        uniq = tuple(sorted(set(seqs)))
        n = len(uniq)
        if n == 0:
            continue
        bset = frozenset(uniq)
        # Small buckets: serial; large: split for worker pool
        if jobs == 1 or n < 5000:
            merged.update(_neighbors_one_length((L, uniq)))
            continue

        n_chunks = min(jobs, max(1, n // 2500))
        chunk_size = (n + n_chunks - 1) // n_chunks
        for i in range(0, n, chunk_size):
            chunk = uniq[i : i + chunk_size]
            chunk_tasks.append((bset, chunk))

    if chunk_tasks:
        with mp.Pool(jobs) as pool:
            for part in pool.imap_unordered(
                _neighbors_chunk_subset, chunk_tasks, chunksize=max(1, len(chunk_tasks) // (jobs * 4))
            ):
                merged.update(part)

    return merged


def redistribute_counts(
    hash_orig: Dict[str, Dict[str, Dict[str, int]]],
    bc1_neighbors: Dict[str, Set[str]],
    bc2_neighbors: Dict[str, Set[str]],
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """Step 3c redistribution: add source counts onto HD1-neighbor pairs that exist in the dict."""
    accum = copy.deepcopy(hash_orig)
    raw_bc_pairs = 0
    alt_bc_pairs = 0
    t0 = time.time()

    for bc1, bc2_dict in hash_orig.items():
        bc1_alts = bc1_neighbors.get(bc1, set())
        for bc2, src in bc2_dict.items():
            raw_bc_pairs += 1
            bc2_alts = bc2_neighbors.get(bc2, set())

            for bc1_alt in bc1_alts:
                if bc2 in accum.get(bc1_alt, {}):
                    alt_bc_pairs += 1
                    for variant, ct in src.items():
                        accum[bc1_alt][bc2][variant] = (
                            accum[bc1_alt][bc2].get(variant, 0) + ct
                        )
                for bc2_alt in bc2_alts:
                    if bc1_alt in accum and bc2_alt in accum[bc1_alt]:
                        alt_bc_pairs += 1
                        for variant, ct in src.items():
                            accum[bc1_alt][bc2_alt][variant] = (
                                accum[bc1_alt][bc2_alt].get(variant, 0) + ct
                            )

            for bc2_alt in bc2_alts:
                if bc1 in accum and bc2_alt in accum[bc1]:
                    alt_bc_pairs += 1
                    for variant, ct in src.items():
                        accum[bc1][bc2_alt][variant] = (
                            accum[bc1][bc2_alt].get(variant, 0) + ct
                        )

            if raw_bc_pairs % 1_000_000 == 0:
                logging.info(
                    "Redistribute progress: raw_pairs=%s alt_ops=%s %.2fs",
                    raw_bc_pairs,
                    alt_bc_pairs,
                    time.time() - t0,
                )
                t0 = time.time()

    logging.info(
        "Redistribution done: raw barcode pairs: %s, alt update ops: %s",
        raw_bc_pairs,
        alt_bc_pairs,
    )
    return accum


def filter_and_write(
    accum: Dict[str, Dict[str, Dict[str, int]]],
    out_path: str,
    min_support_reads: int,
    min_accuracy: float,
) -> None:
    """Write high-confidence pairs as CSV: bc_pair, variant, best_ct, total_ct, accuracy."""
    total_bc_pairs = 0
    valid_bc_pairs = 0
    support_dist: DefaultDict[int, int] = defaultdict(int)

    with open_out(out_path) as fh:
        writer = csv.writer(fh, delimiter=",", lineterminator="\n")
        for bc1 in sorted(accum.keys()):
            for bc2 in sorted(accum[bc1].keys()):
                total_bc_pairs += 1
                var_ct = accum[bc1][bc2]
                best_variant = ""
                best_ct = 0
                total_ct = 0
                for variant in sorted(var_ct.keys()):
                    c = var_ct[variant]
                    total_ct += c
                    if c > best_ct:
                        best_variant = variant
                        best_ct = c

                n_tie = sum(1 for v in var_ct if var_ct[v] == best_ct)
                if n_tie != 1:
                    continue
                support_dist[best_ct] += 1
                acc = best_ct / total_ct if total_ct else 0.0
                if best_ct >= min_support_reads and acc >= min_accuracy:
                    valid_bc_pairs += 1
                    writer.writerow(
                        [
                            f"{bc1}_{bc2}",
                            best_variant,
                            best_ct,
                            total_ct,
                            acc,
                        ]
                    )

    logging.info(
        "Output: total barcode pairs: %s, passing filters: %s",
        total_bc_pairs,
        valid_bc_pairs,
    )
    logging.info("Barcode pair best-support read distribution (unique max):")
    for k in sorted(support_dist.keys()):
        logging.info("  %s\t%s", k, support_dist[k])


def main() -> None:
    args = parse_args()
    jobs = max(1, args.jobs)

    if not os.path.isfile(args.variant_map_in):
        logging.error("Input not found: %s", args.variant_map_in)
        sys.exit(1)

    logging.info("Loading variant-mapped CSV: %s", args.variant_map_in)
    hash_bc, bc1_set, bc2_set, _tl, _vl = load_raw_dict(
        args.variant_map_in,
        ignore_bc1=args.ignore_bc1,
        ignore_bc2=args.ignore_bc2,
    )
    n_pairs = sum(len(v) for v in hash_bc.values())
    logging.info("Raw (bc1, bc2) pairs with data: %s", n_pairs)

    if args.raw_dict_out:
        write_raw_dict_csv(hash_bc, args.raw_dict_out)

    if args.ignore_bc1:
        logging.info("Skipping bc1 HD1 neighbors (--ignore_bc1)")
        nbr1: Dict[str, Set[str]] = {}
    else:
        logging.info("Building HD1 neighbor map for %s bc1 barcodes (%s workers)", len(bc1_set), jobs)
        t0 = time.time()
        nbr1 = build_hd1_neighbors(bc1_set, jobs)
        logging.info("bc1 HD1 map done in %.2fs", time.time() - t0)

    if args.ignore_bc2:
        logging.info("Skipping bc2 HD1 neighbors (--ignore_bc2)")
        nbr2: Dict[str, Set[str]] = {}
    else:
        logging.info("Building HD1 neighbor map for %s bc2 barcodes (%s workers)", len(bc2_set), jobs)
        t0 = time.time()
        nbr2 = build_hd1_neighbors(bc2_set, jobs)
        logging.info("bc2 HD1 map done in %.2fs", time.time() - t0)

    logging.info("Redistributing counts across HD1 neighbors...")
    t0 = time.time()
    accum = redistribute_counts(hash_bc, nbr1, nbr2)
    logging.info("Redistribution finished in %.2fs", time.time() - t0)

    logging.info(
        "Writing filtered dictionary to %s (min_support_reads>=%s, min_accuracy>=%s)",
        args.out,
        args.min_support_reads,
        args.min_accuracy,
    )
    filter_and_write(
        accum,
        args.out,
        args.min_support_reads,
        args.min_accuracy,
    )


if __name__ == "__main__":
    main()
