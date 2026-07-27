#!/usr/bin/env python3
"""
Bar plot of total edit counts per editable template position for one cDNA sample.

Reads a `<rep>.cDNA_parsed.csv` (the step-02 output) and, for each `Pos<N>` column,
sums the per-variant edit counts to get the total number of edit events at that
template position across the whole sample. Bars are drawn in template order; the
primary edit position (from the CSV's `Primary_pos` column, or `--primary_pos` to
override) is drawn in a highlight color.

Usage:
    python plot_edit_positions.py --in_csv <rep>.cDNA_parsed.csv --out_png <path.png>
                                  [--primary_pos <N>] [--title <str>]

The plot is a QC view of *where* editing lands on the template, pooled over variants;
it is not a per-variant figure.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")  # headless: no display on the compute nodes
import matplotlib.pyplot as plt
import pandas as pd


BASE_COLOR = "#4C72B0"       # non-primary positions
PRIMARY_COLOR = "#DD8452"    # primary position (highlighted)


def pos_columns(df: pd.DataFrame):
    """Return [(position:int, column_name:str), ...] in template order for every Pos<N> column."""
    cols = []
    for c in df.columns:
        if c.startswith("Pos"):
            try:
                cols.append((int(c[3:]), c))
            except ValueError:
                continue
    cols.sort(key=lambda t: t[0])
    return cols


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in_csv", required=True, help="Path to a <rep>.cDNA_parsed.csv (step 02 output)")
    parser.add_argument("--out_png", required=True, help="Path for the output PNG")
    parser.add_argument("--primary_pos", type=int, default=None,
                        help="Override the primary position to highlight; default reads it from the CSV")
    parser.add_argument("--title", default=None, help="Plot title; default derived from the CSV filename")
    args = parser.parse_args()

    df = pd.read_csv(args.in_csv)
    cols = pos_columns(df)
    if not cols:
        sys.exit(f"No Pos<N> columns found in {args.in_csv}")

    positions = [p for p, _ in cols]
    totals = [int(df[c].sum()) for _, c in cols]

    # Primary position: explicit flag wins, else the CSV's Primary_pos column (constant per file).
    primary = args.primary_pos
    if primary is None and "Primary_pos" in df.columns and len(df):
        try:
            primary = int(df["Primary_pos"].iloc[0])
        except (ValueError, TypeError):
            primary = None
    if primary is not None and primary not in positions:
        print(f"WARNING: primary position {primary} is not an editable position in this template "
              f"({positions}); no bar will be highlighted", file=sys.stderr)

    colors = [PRIMARY_COLOR if p == primary else BASE_COLOR for p in positions]

    title = args.title or os.path.basename(args.in_csv).replace(".cDNA_parsed.csv", "")

    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(positions) + 2), 4.5))
    x = range(len(positions))
    ax.bar(x, totals, color=colors)
    ax.set_xticks(list(x))
    ax.set_xticklabels([str(p) for p in positions])
    ax.set_xlabel("Editable template position (1-based)")
    ax.set_ylabel("Total edit events (summed over variants)")
    ax.set_title(title)

    if primary in positions:
        # a small legend so the highlight is self-explanatory
        from matplotlib.patches import Patch
        ax.legend(
            handles=[
                Patch(facecolor=PRIMARY_COLOR, label=f"primary position ({primary})"),
                Patch(facecolor=BASE_COLOR, label="other positions"),
            ],
            frameon=False, fontsize=9,
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_png)), exist_ok=True)
    fig.savefig(args.out_png, dpi=150)
    plt.close(fig)
    print(f"edit-position plot -> {args.out_png}")


if __name__ == "__main__":
    main()
