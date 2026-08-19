#!/usr/bin/env python3
"""D-statistic synthetic-mixture calibration curves.

For every outgroup, the synthetic calibration table gives D as a function of an
imposed parent_a (P3) copy number; the observed D of each real test is placed on
that curve to read off a calibrated copy number. One panel per outgroup.
"""
import argparse

import pandas as pd

from plotting import save_figure, write_source_data
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--calibration", required=True, help="synthetic_calibration.tsv")
    p.add_argument("--dstat", required=True, help="dstat.tsv (observed D per test)")
    p.add_argument("--output-base", required=True)
    p.add_argument("--source-data", required=True)
    p.add_argument("--formats", default="png,tiff")
    p.add_argument("--dpi-png", type=int, default=150)
    p.add_argument("--dpi-tiff", type=int, default=300)
    a = p.parse_args()

    cal = pd.read_csv(a.calibration, sep="\t")
    obs = pd.read_csv(a.dstat, sep="\t")
    if cal.empty:
        raise SystemExit(f"no calibration rows in {a.calibration}")
    outgroups = sorted(cal.O.unique())
    n = len(outgroups)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.6), squeeze=False, sharey=True)
    axes = axes[0]
    for ax, out in zip(axes, outgroups):
        c = cal[cal.O == out]
        # calibration curve: use the pooled test if present, else the mean curve
        pooled = c[c.test.str.endswith("_pooled")]
        curve = pooled if not pooled.empty else c.groupby("synthetic_p3_copy_number", as_index=False).D.mean()
        curve = curve.sort_values("synthetic_p3_copy_number")
        ax.plot(curve.synthetic_p3_copy_number, curve.D, color="#333333", lw=1.4, zorder=2,
                label="synthetic calibration")
        o = obs[obs.O == out]
        if not o.empty:
            ax.scatter(o.calibrated_p3_copy_number, o.D, s=26, color="#D55E00",
                       zorder=3, label="observed tests")
            for row in o.itertuples(index=False):
                if pd.notna(row.calibrated_p3_copy_number):
                    ax.vlines(row.calibrated_p3_copy_number, curve.D.min(), row.D,
                              color="#D55E00", lw=0.6, ls=":", zorder=1)
        ax.axhline(0, color="#999999", lw=0.6)
        ax.set_title(f"Outgroup: {out}", loc="left")
        ax.set_xlabel("parent_a (P3) copy number")
        ax.tick_params(length=2)
    axes[0].set_ylabel("D statistic")
    axes[0].legend(loc="best")
    fig.tight_layout()

    formats = [f.strip() for f in a.formats.split(",") if f.strip()]
    save_figure(fig, a.output_base, formats, {"png": a.dpi_png, "tiff": a.dpi_tiff})
    keep = [c for c in ("test", "O", "synthetic_p3_copy_number", "D", "Z_D") if c in cal.columns]
    write_source_data(cal[keep], a.source_data)


if __name__ == "__main__":
    main()
