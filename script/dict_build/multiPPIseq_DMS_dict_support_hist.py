#!/usr/bin/env python3
"""
Distribution of read support per unique barcode.

For each block, counts how many sequencing reads support each unique barcode
(the read-support depth per barcode), then summarizes that distribution. This
quantifies under-sequencing: if most barcodes are seen only once or twice, no
support threshold can yield a large dictionary.

"Support reads" = number of variant-mapped reads carrying that barcode
(pre-redistribution; the raw observed depth). Barcode key honors --ignore_bc2 /
--ignore_bc1 so it matches the dictionary's keying (default --ignore_bc2 = key on
bc1 alone, the intended Dual_NEBuilder design).

Input: step-04 variant_mapped per-read CSVs (read_id = bc1_bc2[_suffix], variant).

Outputs (paths derived from --out_prefix):
  <prefix>.hist.csv      block, support_reads, n_barcodes  (long-format histogram;
                         a final all-blocks "ALL" block is appended)
  <prefix>.summary.csv   per-block + ALL: n_barcodes, n_reads, mean/median/max,
                         and #barcodes (and read fraction) at support >= 1..5,10
  <prefix>.<block>.png   per-block histogram (one combined ALL plot too) if
                         matplotlib is available (skip with --no_plot)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
import sys
from collections import Counter, defaultdict
from typing import Dict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

VM_SUFFIX = "_Dictionary.block_ref_seqs.csv.gz"
THRESHOLDS = (1, 2, 3, 5, 10)


def open_text(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def barcode_support(path: str, ignore_bc1: bool, ignore_bc2: bool) -> Counter:
    """Counter: barcode_key -> number of reads (rows) carrying it."""
    support: Counter = Counter()
    with open_text(path) as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            parts = row[0].split("_")
            if len(parts) < 2:
                continue
            bc1 = "" if ignore_bc1 else parts[0]
            bc2 = "" if ignore_bc2 else parts[1]
            support[(bc1, bc2)] += 1
    return support


def hist_from_support(support: Counter) -> Counter:
    """support-read value -> number of unique barcodes with that value."""
    hist: Counter = Counter()
    for n in support.values():
        hist[n] += 1
    return hist


def summarize(block: str, hist: Counter) -> Dict[str, object]:
    n_barcodes = sum(hist.values())
    n_reads = sum(k * v for k, v in hist.items())
    # median of per-barcode support
    ordered = sorted(hist.items())
    median = 0
    if n_barcodes:
        mid = n_barcodes // 2
        cum = 0
        for k, v in ordered:
            cum += v
            if cum > mid:
                median = k
                break
    row: Dict[str, object] = {
        "block": block,
        "n_unique_barcodes": n_barcodes,
        "n_reads": n_reads,
        "mean_support": round(n_reads / n_barcodes, 3) if n_barcodes else 0,
        "median_support": median,
        "max_support": max(hist) if hist else 0,
    }
    for t in THRESHOLDS:
        bc_ge = sum(v for k, v in hist.items() if k >= t)
        reads_ge = sum(k * v for k, v in hist.items() if k >= t)
        row[f"barcodes_ge{t}"] = bc_ge
        row[f"frac_barcodes_ge{t}"] = (
            round(bc_ge / n_barcodes, 6) if n_barcodes else 0
        )
        row[f"frac_reads_ge{t}"] = round(reads_ge / n_reads, 6) if n_reads else 0
    return row


def plot_hist(block: str, hist: Counter, out_png: str, xmax: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # bin 1..xmax exactly, lump the rest into an ">=xmax+1" bar
    xs = list(range(1, xmax + 1))
    ys = [hist.get(x, 0) for x in xs]
    tail = sum(v for k, v in hist.items() if k > xmax)
    labels = [str(x) for x in xs] + [f"≥{xmax + 1}"]
    ys.append(tail)

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 0.28), 4))
    ax.bar(range(len(ys)), ys, color="#4060c0", width=0.85)
    ax.set_yscale("log")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_xlabel("read support per unique barcode")
    ax.set_ylabel("number of unique barcodes (log)")
    n_bc = sum(hist.values())
    ax.set_title(f"{block}  (n={n_bc:,} barcodes)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def discover_blocks(variant_mapped_dir: str, only):
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
        description="Per-barcode read-support distribution (histogram) -> CSV + PNG."
    )
    p.add_argument(
        "--variant_mapped_dir",
        required=True,
        help="variant_mapped/<LIBTYPE>/blocks (step-04 per-read CSVs)",
    )
    p.add_argument(
        "--out_prefix",
        required=True,
        help="Output path prefix, e.g. .../support_hist (extensions added)",
    )
    p.add_argument("--blocks", nargs="*", default=None, help="Subset of block keys")
    ig = p.add_mutually_exclusive_group()
    ig.add_argument("--ignore_bc1", action="store_true")
    ig.add_argument(
        "--ignore_bc2",
        action="store_true",
        help="Key on bc1 alone (the intended Dual_NEBuilder design)",
    )
    p.add_argument(
        "--xmax",
        type=int,
        default=30,
        help="Plot bins 1..xmax exactly, rest lumped into >=xmax+1 (default 30)",
    )
    p.add_argument("--no_plot", action="store_true", help="Skip PNG plots")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not os.path.isdir(args.variant_mapped_dir):
        logging.error("Not a directory: %s", args.variant_mapped_dir)
        sys.exit(1)
    blocks = discover_blocks(args.variant_mapped_dir, args.blocks)
    if not blocks:
        logging.error("No variant_mapped blocks found.")
        sys.exit(1)

    out_dir = os.path.dirname(os.path.abspath(args.out_prefix))
    os.makedirs(out_dir, exist_ok=True)

    hist_rows = []           # (block, support_reads, n_barcodes)
    summary_rows = []
    combined: Counter = Counter()  # all-blocks histogram

    can_plot = not args.no_plot
    if can_plot:
        try:
            import matplotlib  # noqa: F401
        except Exception as e:  # pragma: no cover
            logging.warning("matplotlib unavailable (%s); skipping plots", e)
            can_plot = False

    for block, vm_path in blocks:
        logging.info("Counting barcode support: %s", block)
        support = barcode_support(vm_path, args.ignore_bc1, args.ignore_bc2)
        hist = hist_from_support(support)
        for k in sorted(hist):
            hist_rows.append((block, k, hist[k]))
            combined[k] += hist[k]
        summary_rows.append(summarize(block, hist))
        if can_plot:
            plot_hist(block, hist, f"{args.out_prefix}.{block}.png", args.xmax)

    # all-blocks combined
    for k in sorted(combined):
        hist_rows.append(("ALL", k, combined[k]))
    summary_rows.append(summarize("ALL", combined))
    if can_plot:
        plot_hist("ALL_blocks", combined, f"{args.out_prefix}.ALL.png", args.xmax)

    hist_csv = f"{args.out_prefix}.hist.csv"
    with open(hist_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block", "support_reads", "n_barcodes"])
        w.writerows(hist_rows)

    summary_csv = f"{args.out_prefix}.summary.csv"
    with open(summary_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    logging.info("Wrote %s and %s", hist_csv, summary_csv)
    if can_plot:
        logging.info("Wrote PNG histograms with prefix %s", args.out_prefix)


if __name__ == "__main__":
    main()
