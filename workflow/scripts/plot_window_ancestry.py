#!/usr/bin/env python3
"""Windowed ancestry: calibrated parent_a copy number along the genome.

Input is ``ancestry/window_ancestry.tsv``. Each target sample gets a panel
showing the calibrated parent_a copy number per fixed-size window, laid out in
genome order with chromosome boundaries marked. Dashed guides mark the integer
copy-number expectations (0..ploidy).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from plotting import natural_key, read_fai_lengths, save_figure, write_source_data
import matplotlib.pyplot as plt


def genome_offsets(chroms, fai_lengths, fallback):
    offset = 0
    offsets = {}
    for c in chroms:
        offsets[c] = offset
        offset += fai_lengths.get(c, fallback.get(c, 0))
    return offsets, offset


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--windows", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--output-base", required=True)
    p.add_argument("--source-data", required=True)
    p.add_argument("--formats", default="png,tiff")
    p.add_argument("--dpi-png", type=int, default=150)
    p.add_argument("--dpi-tiff", type=int, default=300)
    a = p.parse_args()

    w = pd.read_csv(a.windows, sep="\t")
    if w.empty:
        raise SystemExit(f"no windows in {a.windows}")
    fai_lengths = read_fai_lengths(a.fai)
    fai_order = {name: i for i, name in enumerate(fai_lengths)}
    chroms = sorted(w.CHROM.unique(), key=lambda c: (fai_order.get(c, len(fai_order)), natural_key(c)))
    fallback = w.groupby("CHROM").WIN_START.max().to_dict()
    offsets, total = genome_offsets(chroms, fai_lengths, fallback)
    w = w[w.CHROM.isin(offsets)].copy()
    w["genome_mb"] = (w.CHROM.map(offsets) + w.WIN_START) / 1e6

    samples = sorted(w["sample"].unique(), key=natural_key)
    n = len(samples)
    fig, axes = plt.subplots(n, 1, figsize=(11, 1.5 * n), squeeze=False, sharex=True)
    axes = axes[:, 0]
    boundaries = [offsets[c] / 1e6 for c in chroms] + [total / 1e6]
    for ax, sample in zip(axes, samples):
        sd = w[w["sample"] == sample].sort_values("genome_mb")
        ploidy = int(sd.ploidy.iloc[0])
        for k in range(ploidy + 1):
            ax.axhline(k, color="#cccccc", lw=0.6, ls="--", zorder=1)
        for b in boundaries:
            ax.axvline(b, color="#eeeeee", lw=0.6, zorder=0)
        ax.plot(sd.genome_mb, sd.calibrated_parent_a_copy_number, color="#0072B2",
                lw=0.9, marker="o", ms=2, zorder=3)
        ax.set_ylim(-0.2, ploidy + 0.2)
        ax.set_yticks(range(ploidy + 1))
        ax.set_ylabel("parent_a copies")
        ax.set_title(sample, loc="left", fontweight="bold")
    tick_centers = [(boundaries[i] + boundaries[i + 1]) / 2 for i in range(len(chroms))]
    axes[-1].set_xticks(tick_centers)
    axes[-1].set_xticklabels(chroms, rotation=90)
    axes[-1].set_xlabel("Chromosome (genome order)")
    fig.tight_layout()

    formats = [f.strip() for f in a.formats.split(",") if f.strip()]
    save_figure(fig, a.output_base, formats, {"png": a.dpi_png, "tiff": a.dpi_tiff})
    write_source_data(
        w[["sample", "CHROM", "WIN_START", "n_sites", "calibrated_parent_a_fraction",
           "calibrated_parent_a_copy_number"]],
        a.source_data,
    )


if __name__ == "__main__":
    main()
