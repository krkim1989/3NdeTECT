#!/usr/bin/env python3
"""nQuire ploidy evidence: allele-spectrum thirds-to-half ratio per sample.

Input is ``ploidy/ploidy_summary.tsv``. Each sample is a bar of its
thirds-to-half allele-balance density ratio (triploid signal / diploid signal);
the dashed line is the diploid-control-derived decision threshold for that
mapping panel, and bars are coloured by the resulting ploidy call.
"""
import argparse

import pandas as pd

from plotting import natural_key, save_figure, write_source_data
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

CALL_COLORS = {
    "triploid_supported": "#D55E00",
    "triploid_probable": "#E69F00",
    "triploid_not_supported": "#0072B2",
    "nontriploid_control": "#56B4E9",
    "insufficient_data": "#BBBBBB",
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ploidy", required=True, help="ploidy_summary.tsv")
    p.add_argument("--output-base", required=True)
    p.add_argument("--source-data", required=True)
    p.add_argument("--formats", default="png,tiff")
    p.add_argument("--dpi-png", type=int, default=150)
    p.add_argument("--dpi-tiff", type=int, default=300)
    a = p.parse_args()

    d = pd.read_csv(a.ploidy, sep="\t")
    if d.empty:
        raise SystemExit(f"no ploidy rows in {a.ploidy}")
    order = sorted(d["sample"].unique(), key=natural_key)
    d = d.set_index("sample").loc[order].reset_index()
    colors = [CALL_COLORS.get(c, "#777777") for c in d.ploidy_call]

    x = range(len(d))
    fig, ax = plt.subplots(figsize=(max(5, 0.55 * len(d) + 2), 4))
    ax.bar(x, d.thirds_to_half_ratio.to_numpy(float), color=colors, width=0.75,
           edgecolor="white", linewidth=0.3)
    # one threshold line per distinct (panel, threshold)
    for thr in sorted(d.diploid_ratio_threshold.dropna().unique()):
        ax.axhline(thr, color="#333333", lw=0.9, ls="--")
        ax.text(len(d) - 0.5, thr, f" threshold={thr:.3g}", va="bottom", ha="right", fontsize=7)
    ax.set_xticks(list(x))
    ax.set_xticklabels(d["sample"], rotation=90)
    ax.set_ylabel("thirds-to-half density ratio")
    ax.set_xlabel("Sample")
    present = [c for c in CALL_COLORS if c in set(d.ploidy_call)]
    handles = [mpatches.Patch(color=CALL_COLORS[c], label=c.replace("_", " ")) for c in present]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()

    formats = [f.strip() for f in a.formats.split(",") if f.strip()]
    save_figure(fig, a.output_base, formats, {"png": a.dpi_png, "tiff": a.dpi_tiff})
    keep = ["sample", "group", "role", "mapping_panel", "expected_ploidy", "best_model",
            "thirds_to_half_ratio", "diploid_ratio_threshold", "ploidy_call"]
    write_source_data(d[[c for c in keep if c in d.columns]], a.source_data)


if __name__ == "__main__":
    main()
