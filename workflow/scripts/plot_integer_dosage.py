#!/usr/bin/env python3
"""Integer dosage-state composition per target sample.

Input is ``hybrid_validation/dosage_states.tsv``. For each target the fraction
of diagnostic sites assigned to each integer parent_a copy state (0..ploidy) is
drawn as a stacked bar, so a clean triploid backcross shows as a dominant single
state while mosaic/aneuploid patterns spread across states.
"""
import argparse

import pandas as pd

from plotting import copy_color, copy_label, natural_key, save_figure, write_source_data
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--states", required=True, help="dosage_states.tsv")
    p.add_argument("--ancestry-label", default="parent_a")
    p.add_argument("--output-base", required=True)
    p.add_argument("--source-data", required=True)
    p.add_argument("--formats", default="png,tiff")
    p.add_argument("--dpi-png", type=int, default=150)
    p.add_argument("--dpi-tiff", type=int, default=300)
    a = p.parse_args()

    d = pd.read_csv(a.states, sep="\t")
    if d.empty:
        raise SystemExit(f"no dosage-state rows in {a.states}")
    order = sorted(d["sample"].unique(), key=natural_key)
    d = d.set_index("sample").loc[order].reset_index()
    max_ploidy = int(d.ploidy.max())
    frac_cols = [f"state_{k}_fraction" for k in range(max_ploidy + 1)]
    for col in frac_cols:
        if col not in d.columns:
            d[col] = 0.0
    d[frac_cols] = d[frac_cols].fillna(0.0)

    samples = d["sample"].tolist()
    x = range(len(samples))
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(samples) + 2), 4))
    bottom = [0.0] * len(samples)
    for k in range(max_ploidy + 1):
        vals = d[f"state_{k}_fraction"].to_numpy(float)
        ax.bar(x, vals, bottom=bottom, color=copy_color(k), width=0.75,
               edgecolor="white", linewidth=0.3)
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_xticks(list(x))
    ax.set_xticklabels(samples, rotation=90)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of diagnostic sites")
    ax.set_xlabel("Target sample")
    handles = [mpatches.Patch(color=copy_color(k), label=copy_label(k, a.ancestry_label))
               for k in range(max_ploidy + 1)]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()

    formats = [f.strip() for f in a.formats.split(",") if f.strip()]
    save_figure(fig, a.output_base, formats, {"png": a.dpi_png, "tiff": a.dpi_tiff})
    write_source_data(d[["sample", "ploidy", "n_sites", "pct_near_integer_dosage"] + frac_cols],
                      a.source_data)


if __name__ == "__main__":
    main()
