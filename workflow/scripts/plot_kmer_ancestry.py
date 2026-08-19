#!/usr/bin/env python3
"""Mapping-independent k-mer ancestry copy number per sample.

Input is ``kmer/kmer_ancestry.tsv``. Bars show the calibrated parent_a copy
number estimated purely from diagnostic k-mer counts, coloured by sample role,
providing an assembly/mapping-independent cross-check of the read-mapped dosage.
"""
import argparse

import pandas as pd

from plotting import natural_key, save_figure, write_source_data
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

ROLE_COLORS = {
    "target": "#D55E00",
    "parent_a": "#0072B2",
    "parent_b": "#E69F00",
    "outgroup": "#999999",
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kmer", required=True, help="kmer_ancestry.tsv")
    p.add_argument("--output-base", required=True)
    p.add_argument("--source-data", required=True)
    p.add_argument("--formats", default="png,tiff")
    p.add_argument("--dpi-png", type=int, default=150)
    p.add_argument("--dpi-tiff", type=int, default=300)
    a = p.parse_args()

    d = pd.read_csv(a.kmer, sep="\t")
    if d.empty:
        raise SystemExit(f"no k-mer rows in {a.kmer}")
    order = sorted(d["sample"].unique(), key=natural_key)
    d = d.set_index("sample").loc[order].reset_index()
    colors = [ROLE_COLORS.get(r, "#777777") for r in d.role]

    x = range(len(d))
    fig, ax = plt.subplots(figsize=(max(5, 0.55 * len(d) + 2), 4))
    ax.bar(x, d.calibrated_parent_a_copy_number.to_numpy(float), color=colors, width=0.75,
           edgecolor="white", linewidth=0.3)
    ax.set_xticks(list(x))
    ax.set_xticklabels(d["sample"], rotation=90)
    ax.set_ylabel("k-mer calibrated parent_a copies")
    ax.set_xlabel("Sample")
    present = [r for r in ROLE_COLORS if r in set(d.role)]
    handles = [mpatches.Patch(color=ROLE_COLORS[r], label=r) for r in present]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()

    formats = [f.strip() for f in a.formats.split(",") if f.strip()]
    save_figure(fig, a.output_base, formats, {"png": a.dpi_png, "tiff": a.dpi_tiff})
    keep = ["sample", "group", "role", "ploidy", "raw_parent_a_fraction",
            "calibrated_parent_a_fraction", "calibrated_parent_a_copy_number"]
    write_source_data(d[[c for c in keep if c in d.columns]], a.source_data)


if __name__ == "__main__":
    main()
