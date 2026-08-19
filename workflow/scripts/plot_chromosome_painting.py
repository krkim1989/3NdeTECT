#!/usr/bin/env python3
"""Chromosome painting: per-target ancestry copy number along each chromosome.

Input is the local-ancestry HMM segment table
(``ancestry/local_ancestry_segments.tsv``). Each target sample becomes a
stacked panel; within a panel every chromosome is a horizontal track drawn to
its true length from the reference .fai, with a grey background for stretches
that carry no diagnostic markers and coloured tracts where the HMM assigned an
integer parent_a copy number.

Fully generic: the number of chromosomes, samples, and ploidy states are taken
from the data, and the ancestry name shown in the legend comes from the config.
"""
import argparse
from pathlib import Path

import pandas as pd

from plotting import (
    NO_DATA_COLOR,
    copy_color,
    copy_label,
    natural_key,
    read_fai_lengths,
    save_figure,
    write_source_data,
)
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


def ordered_chromosomes(chroms, fai_lengths):
    """Chromosomes present in the data, in reference (.fai) then natural order."""
    fai_order = {name: i for i, name in enumerate(fai_lengths)}
    return sorted(chroms, key=lambda c: (fai_order.get(c, len(fai_order)), natural_key(c)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--segments", required=True, help="local_ancestry_segments.tsv")
    p.add_argument("--fai", required=True, help="reference .fai for chromosome lengths")
    p.add_argument("--ancestry-label", default="parent_a", help="legend name for the copied ancestry")
    p.add_argument("--focal-copy", choices=["parent_a", "parent_b"], default="parent_a",
                   help="which parent's copy number to paint; parent_b uses ploidy - A_COPY")
    p.add_argument("--model", default="constrained", help="HMM model column to paint")
    p.add_argument("--output-base", required=True, help="output path without extension")
    p.add_argument("--source-data", required=True)
    p.add_argument("--formats", default="png,tiff")
    p.add_argument("--dpi-png", type=int, default=150)
    p.add_argument("--dpi-tiff", type=int, default=300)
    a = p.parse_args()

    seg = pd.read_csv(a.segments, sep="\t")
    if "model" in seg.columns:
        seg = seg[seg.model == a.model].copy()
    if seg.empty:
        raise SystemExit(f"no '{a.model}' segments to paint in {a.segments}")

    # Paint from the perspective of the requested parent. parent_b copies are the
    # complement of parent_a copies within each site's ploidy, so both species can
    # be painted from the same segment table.
    if a.focal_copy == "parent_b":
        if "ploidy" not in seg.columns:
            raise SystemExit("parent_b painting needs a 'ploidy' column in the segments table")
        seg["PAINT_COPY"] = (seg.ploidy - seg.A_COPY).astype(int)
    else:
        seg["PAINT_COPY"] = seg.A_COPY.astype(int)

    fai_lengths = read_fai_lengths(a.fai)
    samples = sorted(seg["sample"].unique(), key=natural_key)
    chroms = ordered_chromosomes(seg.CHROM.unique(), fai_lengths)
    chrom_index = {c: i for i, c in enumerate(chroms)}
    max_copy = int(seg.PAINT_COPY.max())

    # length fallback: if a contig is missing from the .fai, use the furthest
    # segment end so the track is still drawn to a sensible extent.
    seg_max_end = seg.groupby("CHROM").end.max().to_dict()
    lengths = {c: fai_lengths.get(c, seg_max_end.get(c, 1)) for c in chroms}
    genome_max_mb = max(lengths.values()) / 1e6

    n = len(samples)
    panel_h = 0.28 * len(chroms) + 0.5  # inches per sample panel, scales with chromosome count
    fig, axes = plt.subplots(
        n, 1, figsize=(11, panel_h * n), squeeze=False, sharex=True,
    )
    axes = axes[:, 0]
    bar_h = 0.7
    for ax, sample in zip(axes, samples):
        sd = seg[seg["sample"] == sample]
        for chrom in chroms:
            y = chrom_index[chrom]
            ax.broken_barh([(0, lengths[chrom] / 1e6)], (y - bar_h / 2, bar_h),
                           facecolors=NO_DATA_COLOR, edgecolors="none", zorder=1)
            for row in sd[sd.CHROM == chrom].itertuples(index=False):
                x0 = row.start / 1e6
                width = max((row.end - row.start) / 1e6, genome_max_mb * 0.0015)
                ax.broken_barh([(x0, width)], (y - bar_h / 2, bar_h),
                               facecolors=copy_color(row.PAINT_COPY), edgecolors="none", zorder=2)
        ax.set_yticks(range(len(chroms)))
        ax.set_yticklabels(chroms)
        ax.set_ylim(-0.6, len(chroms) - 0.4)
        ax.set_xlim(0, genome_max_mb * 1.01)
        ax.set_ylabel("Chromosome")
        ax.set_title(sample, loc="left", fontweight="bold")
        ax.tick_params(length=2)
    axes[-1].set_xlabel("Position (Mbp)")

    handles = [mpatches.Patch(color=copy_color(k), label=copy_label(k, a.ancestry_label))
               for k in range(max_copy + 1)]
    handles.append(mpatches.Patch(color=NO_DATA_COLOR, label="no diagnostic markers"))
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 5),
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    formats = [f.strip() for f in a.formats.split(",") if f.strip()]
    save_figure(fig, a.output_base, formats, {"png": a.dpi_png, "tiff": a.dpi_tiff})
    src_cols = [c for c in ["sample", "CHROM", "start", "end", "ploidy", "A_COPY", "PAINT_COPY", "n_sites"]
                if c in seg.columns]
    write_source_data(seg[src_cols], a.source_data)


if __name__ == "__main__":
    main()
