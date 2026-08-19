import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"


def load(name):
    p = SCRIPTS / name
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_copy_label_singular_and_plural():
    m = load("plotting.py")
    assert m.copy_label(1, "C. gigas") == "1 C. gigas copy"
    assert m.copy_label(2, "C. gigas") == "2 C. gigas copies"


def test_natural_key_orders_chromosomes_numerically():
    m = load("plotting.py")
    names = ["chr10", "chr2", "chr1"]
    assert sorted(names, key=m.natural_key) == ["chr1", "chr2", "chr10"]


def test_save_figure_writes_every_requested_format(tmp_path):
    m = load("plotting.py")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    written = m.save_figure(fig, tmp_path / "demo", ["png", "tiff"], {"png": 90, "tiff": 120})
    assert {p.suffix for p in written} == {".png", ".tiff"}
    assert all(p.is_file() and p.stat().st_size > 0 for p in written)


def _fai(tmp_path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t3000000\t0\t60\t61\nchr2\t2000000\t3000000\t60\t61\n")
    return fai


def test_chromosome_painting_renders_png_tiff_and_source(tmp_path):
    seg = pd.DataFrame([
        dict(sample="S1", model="constrained", CHROM="chr1", start=0, end=1500000,
             bp_length=1500000, A_COPY=1, n_sites=10, mean_A_read_fraction=0.33, pass_filter=True),
        dict(sample="S1", model="constrained", CHROM="chr1", start=1500000, end=3000000,
             bp_length=1500000, A_COPY=2, n_sites=8, mean_A_read_fraction=0.66, pass_filter=True),
        dict(sample="S1", model="constrained", CHROM="chr2", start=0, end=2000000,
             bp_length=2000000, A_COPY=0, n_sites=12, mean_A_read_fraction=0.0, pass_filter=True),
    ])
    seg_path = tmp_path / "segments.tsv"
    seg.to_csv(seg_path, sep="\t", index=False)
    base = tmp_path / "chromosome_painting"
    src = tmp_path / "source.tsv"
    subprocess.run([
        sys.executable, str(SCRIPTS / "plot_chromosome_painting.py"),
        "--segments", str(seg_path), "--fai", str(_fai(tmp_path)),
        "--ancestry-label", "C. gigas", "--output-base", str(base),
        "--source-data", str(src), "--formats", "png,tiff",
    ], cwd=ROOT, check=True, capture_output=True, text=True)
    assert base.with_suffix(".png").is_file()
    assert base.with_suffix(".tiff").is_file()
    assert not pd.read_csv(src, sep="\t").empty


def test_parent_b_painting_uses_complement_copy_number(tmp_path):
    # ploidy 3, parent_a copies 1 -> parent_b copies must be painted as 2
    seg = pd.DataFrame([
        dict(sample="S1", model="constrained", CHROM="chr1", start=0, end=3000000,
             bp_length=3000000, ploidy=3, A_COPY=1, n_sites=10,
             mean_A_read_fraction=0.33, pass_filter=True),
    ])
    seg_path = tmp_path / "segments.tsv"
    seg.to_csv(seg_path, sep="\t", index=False)
    base = tmp_path / "chromosome_painting_parentB"
    src = tmp_path / "source_b.tsv"
    subprocess.run([
        sys.executable, str(SCRIPTS / "plot_chromosome_painting.py"),
        "--segments", str(seg_path), "--fai", str(_fai(tmp_path)),
        "--ancestry-label", "species B", "--focal-copy", "parent_b",
        "--output-base", str(base), "--source-data", str(src), "--formats", "png",
    ], cwd=ROOT, check=True, capture_output=True, text=True)
    assert base.with_suffix(".png").is_file()
    out = pd.read_csv(src, sep="\t")
    assert int(out.PAINT_COPY.iloc[0]) == 2


def _run_painting(tmp_path, seg, fai_text, extra=()):
    seg_path = tmp_path / "segments.tsv"
    seg.to_csv(seg_path, sep="\t", index=False)
    fai = tmp_path / "ref.fa.fai"
    fai.write_text(fai_text)
    base = tmp_path / "painting"
    src = tmp_path / "source.tsv"
    return subprocess.run([
        sys.executable, str(SCRIPTS / "plot_chromosome_painting.py"),
        "--segments", str(seg_path), "--fai", str(fai),
        "--output-base", str(base), "--source-data", str(src), "--formats", "png", *extra,
    ], cwd=ROOT, check=True, capture_output=True, text=True)


def test_contig_name_mismatch_warns(tmp_path):
    # a segment chromosome absent from the .fai must trigger a stderr warning
    seg = pd.DataFrame([
        dict(sample="S1", model="constrained", CHROM="scaffold_x", start=0, end=500000,
             bp_length=500000, ploidy=3, A_COPY=0, n_sites=5, mean_A_read_fraction=0.0,
             pass_filter=True),
    ])
    r = _run_painting(tmp_path, seg, "chr1\t10000000\t0\t60\t61\n")
    assert "not found in the .fai" in r.stderr


def test_all_fai_chromosomes_includes_markerless_contigs(tmp_path):
    # chr2 has no segments but must still be painted in --all-fai-chromosomes mode
    seg = pd.DataFrame([
        dict(sample="S1", model="constrained", CHROM="chr1", start=0, end=3000000,
             bp_length=3000000, ploidy=3, A_COPY=1, n_sites=10, mean_A_read_fraction=0.33,
             pass_filter=True),
    ])
    fai = "chr1\t10000000\t0\t60\t61\nchr2\t8000000\t10000000\t60\t61\n"
    r = _run_painting(tmp_path, seg, fai, extra=["--all-fai-chromosomes"])
    assert (tmp_path / "painting.png").is_file()
    assert r.returncode == 0
