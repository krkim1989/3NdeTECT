#!/usr/bin/env python3
"""Shared plotting helpers for the 3NdeTECT figures module.

Every figure script imports from here so that the whole figure set reads as one
system: identical fonts, a single copy-number colour scheme, and a uniform
"save PNG + publication TIFF + machine-readable source-data TSV" contract.

Nothing here is organism-specific. Ancestry labels (e.g. the parent_a group
name) are always passed in from the caller, which reads them from the run
config, so the same code paints an oyster, a loach, or any other cross.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ---------------------------------------------------------------------------
# House style
# ---------------------------------------------------------------------------
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "figure.dpi": 100,
        "savefig.bbox": "tight",
    }
)

# Okabe-Ito colour-blind-safe palette, ordered so that copy number 0,1,2,...
# get distinct, high-contrast hues. Index by integer copy number; anything
# beyond the list cycles. Missing/uncovered regions use NO_DATA_COLOR.
COPY_PALETTE = [
    "#0072B2",  # 0 copies - blue
    "#E69F00",  # 1 copy   - orange
    "#D55E00",  # 2 copies - vermillion
    "#009E73",  # 3 copies - bluish green
    "#CC79A7",  # 4 copies - reddish purple
    "#56B4E9",  # 5 copies - sky blue
    "#F0E442",  # 6 copies - yellow
]
NO_DATA_COLOR = "#DDDDDD"


def copy_color(copy_number: int) -> str:
    """Discrete colour for an integer ancestry copy number."""
    return COPY_PALETTE[int(copy_number) % len(COPY_PALETTE)]


def copy_label(copy_number: int, ancestry_label: str) -> str:
    """Legend text, e.g. '1 C. gigas copy' / '2 C. gigas copies'."""
    unit = "copy" if int(copy_number) == 1 else "copies"
    return f"{int(copy_number)} {ancestry_label} {unit}"


def natural_key(name: str):
    """Sort key that orders chr1, chr2, ..., chr10 the way humans expect."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(name))]


def read_fai_lengths(fai_path: str | Path) -> dict[str, int]:
    """Return {contig: length} from a samtools .fai index."""
    lengths: dict[str, int] = {}
    for line in Path(fai_path).read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        lengths[fields[0]] = int(fields[1])
    return lengths


def save_figure(fig, out_base: str | Path, formats, dpi) -> list[Path]:
    """Save ``fig`` to every requested raster format next to ``out_base``.

    ``out_base`` is a path without extension. ``formats`` is a list such as
    ["png", "tiff"]. ``dpi`` is either a single int or a {format: dpi} mapping,
    so PNGs can stay light (screen/HTML) while TIFFs are print resolution.
    TIFFs are LZW-compressed (needs Pillow) to keep publication files small.
    """
    out_base = Path(out_base)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for fmt in formats:
        resolution = dpi.get(fmt, 300) if isinstance(dpi, dict) else int(dpi)
        target = out_base.with_suffix(f".{fmt}")
        tmp = target.with_suffix(f".{fmt}.tmp")
        save_kwargs = {"dpi": resolution, "format": fmt}
        if fmt in ("tiff", "tif"):
            save_kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(tmp, **save_kwargs)
        tmp.replace(target)
        written.append(target)
    plt.close(fig)
    return written


def write_source_data(df, path: str | Path) -> None:
    """Persist the exact numbers behind a figure as a TSV (atomic write)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format="%.7g")
    tmp.replace(path)
