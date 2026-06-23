#!/usr/bin/env python3
"""
Per-barcode support x fidelity scatter, to visually assess a dictionary.

For each block we load + redistribute the barcode counts ONCE (reusing the step3
functions, so per-barcode numbers are identical to the real pipeline), then for
every barcode (pair) derive:
  best_ct   = reads supporting its most-supported variant (post-redistribution)
  total_ct  = all reads on that barcode (post-redistribution)
  accuracy  = best_ct / total_ct          (a.k.a. fidelity)

It then plots one dot per barcode on log2-log2 axes:
  x = best_ct (support reads)   ticks 2,4,8,...,512
  y = accuracy (fidelity)       ticks 1/64, 1/32, ..., 1/2, 1
Dots are colored by local barcode density (overplotting is severe; a raw scatter
of ~200k identical grid points would be unreadable).

This is an EXPLORATORY diagnostic: you read the cloud to CHOOSE per-block support/
fidelity cutoffs. It is therefore produced independently of the decided parameters
-- the dashed reference lines (--min_support/--min_accuracy, default 3 x 0.9) are a
generic guide for reading the plot, NOT the finally-chosen cutoffs.

Tie barcodes (no unique best variant) are kept for the cloud but flagged via the
n_tie column in the per-barcode CSV; pass --drop_ties to exclude them (matching
the dictionary filter, which never calls a tied barcode).

--ignore_bc1/--ignore_bc2 must match how step 05 was run (default: --ignore_bc2)
so the cloud reflects the actual dictionary.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import logging
import os
import sys
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

VM_SUFFIX = "_Dictionary.block_ref_seqs.csv.gz"

# Axis tick targets requested for the assessment plot.
X_TICKS = [2, 4, 8, 16, 32, 64, 128, 256, 512]
Y_TICKS = [1 / 64, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0]

PERBC_FIELDS = ["block", "bc_pair", "best_ct", "total_ct", "accuracy", "n_tie"]


def load_step3(script_dir: str):
    path = os.path.join(script_dir, "multiPPIseq_DMS_dict_build_step3.py")
    spec = importlib.util.spec_from_file_location("step3", path)
    mod = importlib.util.module_from_spec(spec)
    # Register before exec so step3's multiprocessing workers (pickled by
    # qualified name "step3.<fn>") resolve in worker processes.
    sys.modules["step3"] = mod
    spec.loader.exec_module(mod)
    return mod


def per_barcode(
    step3,
    variant_map_in: str,
    *,
    ignore_bc1: bool,
    ignore_bc2: bool,
    jobs: int,
) -> List[Tuple[str, int, int, float, int]]:
    """One tuple per barcode pair: (bc_pair, best_ct, total_ct, accuracy, n_tie)."""
    hash_bc, bc1_set, bc2_set, _t, _v = step3.load_raw_dict(
        variant_map_in, ignore_bc1=ignore_bc1, ignore_bc2=ignore_bc2
    )
    nbr1 = {} if ignore_bc1 else step3.build_hd1_neighbors(bc1_set, jobs)
    nbr2 = {} if ignore_bc2 else step3.build_hd1_neighbors(bc2_set, jobs)
    accum = step3.redistribute_counts(hash_bc, nbr1, nbr2)

    out: List[Tuple[str, int, int, float, int]] = []
    for bc1, bc2_dict in accum.items():
        for bc2, var_ct in bc2_dict.items():
            best_ct = 0
            total_ct = 0
            for c in var_ct.values():
                total_ct += c
                if c > best_ct:
                    best_ct = c
            n_tie = sum(1 for c in var_ct.values() if c == best_ct)
            acc = best_ct / total_ct if total_ct else 0.0
            out.append((f"{bc1}_{bc2}", best_ct, total_ct, acc, n_tie))
    return out


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


def density_color(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Local barcode count per point via a 2D histogram of log2 coords."""
    lx = np.log2(x)
    ly = np.log2(y)
    nb = 80
    H, xe, ye = np.histogram2d(lx, ly, bins=nb)
    ix = np.clip(np.digitize(lx, xe) - 1, 0, nb - 1)
    iy = np.clip(np.digitize(ly, ye) - 1, 0, nb - 1)
    return H[ix, iy]


def plot_scatter(
    title: str,
    best_ct: np.ndarray,
    acc: np.ndarray,
    out_png: str,
    *,
    min_support: int,
    min_accuracy: float,
) -> None:
    if best_ct.size == 0:
        logging.warning("No barcodes for %s; skipping plot", title)
        return

    # accuracy floor for the log axis: lowest requested tick, or data min if lower.
    y_floor = min(Y_TICKS[0], float(acc[acc > 0].min()) if np.any(acc > 0) else Y_TICKS[0])
    y = np.clip(acc, y_floor, 1.0)
    x = np.clip(best_ct, 1, None).astype(float)

    dens = density_color(x, y)
    order = np.argsort(dens)  # densest last (drawn on top)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    sc = ax.scatter(
        x[order],
        y[order],
        c=dens[order],
        s=5,
        alpha=0.35,
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(2.0, dens.max())),
        linewidths=0,
        rasterized=True,
    )
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("barcode density (count in 2D bin)")

    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(X_TICKS)
    ax.set_xticklabels([str(t) for t in X_TICKS])
    ax.set_yticks(Y_TICKS)
    ax.set_yticklabels(["1/64", "1/32", "1/16", "1/8", "1/4", "1/2", "1"])
    ax.set_xlim(1, max(X_TICKS[-1], float(x.max()) * 1.1))
    ax.set_ylim(y_floor * 0.95, 1.05)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

    # reference threshold lines (a generic guide for reading the cloud — NOT the
    # finally-chosen per-block cutoffs; those are decided by inspecting this plot)
    ax.axvline(min_support, color="crimson", ls="--", lw=1)
    ax.axhline(min_accuracy, color="crimson", ls="--", lw=1)
    ax.axvspan(min_support, ax.get_xlim()[1], ymin=0, ymax=1, color="crimson", alpha=0.04)
    ax.set_title(f"{title} — {best_ct.size:,} barcodes")
    ax.set_xlabel("support reads for best variant (log2)")
    ax.set_ylabel("fidelity = best_ct / total_ct (log2)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    logging.info("Wrote %s", out_png)


def parse_args():
    p = argparse.ArgumentParser(
        description="Per-barcode support x fidelity scatter for dictionary assessment."
    )
    p.add_argument("--variant_mapped_dir", required=True)
    p.add_argument("--out_prefix", required=True, help="Output path prefix")
    p.add_argument(
        "--script_dir", default=os.path.dirname(os.path.abspath(__file__))
    )
    p.add_argument("--blocks", nargs="*", default=None)
    p.add_argument(
        "--min_support",
        type=int,
        default=3,
        help="Reference line x (default 3). A generic guide only, not a chosen cutoff.",
    )
    p.add_argument(
        "--min_accuracy",
        type=float,
        default=0.9,
        help="Reference line y (default 0.9). A generic guide only, not a chosen cutoff.",
    )
    p.add_argument("--drop_ties", action="store_true", help="Exclude tie barcodes (match dict filter)")
    ig = p.add_mutually_exclusive_group()
    ig.add_argument("--ignore_bc1", action="store_true")
    ig.add_argument("--ignore_bc2", action="store_true")
    p.add_argument("--jobs", "-j", type=int, default=1)
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

    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)) or ".", exist_ok=True)
    csv_path = f"{args.out_prefix}.perbarcode.csv.gz"
    fh = gzip.open(csv_path, "wt", newline="")
    writer = csv.writer(fh, lineterminator="\n")
    writer.writerow(PERBC_FIELDS)

    all_best: List[int] = []
    all_acc: List[float] = []
    for block, vm_path in blocks:
        logging.info("Block %s: loading + redistributing once", block)
        rows = per_barcode(
            step3, vm_path, ignore_bc1=args.ignore_bc1, ignore_bc2=args.ignore_bc2, jobs=max(1, args.jobs)
        )
        b_best, b_acc = [], []
        for bc_pair, best_ct, total_ct, acc, n_tie in rows:
            writer.writerow([block, bc_pair, best_ct, total_ct, f"{acc:.6g}", n_tie])
            if args.drop_ties and n_tie != 1:
                continue
            b_best.append(best_ct)
            b_acc.append(acc)
        logging.info("Block %s: %s barcodes plotted", block, len(b_best))
        plot_scatter(
            block,
            np.asarray(b_best),
            np.asarray(b_acc),
            f"{args.out_prefix}.{block}.png",
            min_support=args.min_support,
            min_accuracy=args.min_accuracy,
        )
        all_best.extend(b_best)
        all_acc.extend(b_acc)

    fh.close()
    logging.info("Per-barcode table -> %s", csv_path)

    plot_scatter(
        "ALL blocks",
        np.asarray(all_best),
        np.asarray(all_acc),
        f"{args.out_prefix}.ALL.png",
        min_support=args.min_support,
        min_accuracy=args.min_accuracy,
    )


if __name__ == "__main__":
    main()
