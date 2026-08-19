#!/usr/bin/env python3
import argparse,base64,html
from pathlib import Path
import pandas as pd

STYLE="""body{font-family:system-ui,sans-serif;max-width:1200px;margin:2rem auto;color:#183047}h1,h2{color:#173b5e}table{border-collapse:collapse;width:100%;margin-bottom:2rem}th,td{border:1px solid #ccd6df;padding:.45rem;text-align:right}th{background:#e7f0f7}td:first-child,th:first-child{text-align:left}.warn{background:#fff2cc;padding:1rem;border-left:5px solid #e69f00}.ok{background:#e2f1e8;padding:1rem;border-left:5px solid #2c7e57}pre{background:#f1f4f6;padding:1rem;overflow:auto}figure{margin:0 0 2rem}figure img{max-width:100%;height:auto;border:1px solid #dde6ee}figcaption{font-size:.85em;color:#4a5b6b;margin-top:.3rem}"""

# Figure stem -> caption, in report display order. Missing files are skipped, so
# optional figures (e.g. k-mer) simply do not appear.
FIGURE_CAPTIONS=[
    ("chromosome_painting","Chromosome painting: per-target parent_a copy number along each chromosome."),
    ("chromosome_painting_parentB","Chromosome painting from the parent_b perspective (ploidy - parent_a copies)."),
    ("window_ancestry","Windowed ancestry: calibrated parent_a copy number in genome order."),
    ("integer_dosage","Integer dosage-state composition per target."),
    ("dstat_calibration","D-statistic synthetic-mixture calibration and observed tests."),
    ("nquire_ploidy","nQuire allele-spectrum ploidy evidence per sample."),
    ("kmer_ancestry","Mapping-independent k-mer ancestry copy number."),
]

def table(path):
    d=pd.read_csv(path,sep="\t")
    if d.empty: raise ValueError(f"empty report input: {path}")
    return d.to_html(index=False,float_format=lambda x:f"{x:.4g}",escape=True)

def figures_section(figures_dir):
    if not figures_dir: return ""
    d=Path(figures_dir); blocks=[]
    for stem,caption in FIGURE_CAPTIONS:
        png=d/f"{stem}.png"
        if not png.is_file(): continue
        uri=base64.b64encode(png.read_bytes()).decode("ascii")
        blocks.append(f"<figure><img alt='{html.escape(stem)}' src='data:image/png;base64,{uri}'>"
                      f"<figcaption>{html.escape(caption)}</figcaption></figure>")
    return f"<h2>Figures</h2>{''.join(blocks)}" if blocks else ""

def main():
    p=argparse.ArgumentParser();
    for name in ("hybrid","ploidy","dosage","states","local","windows","dstat","dcal","config","versions"): p.add_argument(f"--{name}",required=True)
    p.add_argument("--reciprocal"); p.add_argument("--kmer"); p.add_argument("--cox1"); p.add_argument("--figures-dir")
    p.add_argument("--title",default="3NdeTECT report"); p.add_argument("--output",required=True); a=p.parse_args()
    cfg=html.escape(Path(a.config).read_text())
    figures=figures_section(a.figures_dir)
    reciprocal=f"<h2>Reciprocal-reference concordance</h2>{table(a.reciprocal)}" if a.reciprocal else ""
    kmer=f"<h2>Mapping-independent k-mer ancestry</h2>{table(a.kmer)}" if a.kmer else ""
    cox1=f"<h2>FASTQ-derived COX1 consensus and maternal lineage</h2>{table(a.cox1)}" if a.cox1 else ""
    body=f"""<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(a.title)}</title><style>{STYLE}</style></head><body><h1>{html.escape(a.title)}</h1><div class='warn'><b>Interpretation:</b> ancestry dosage and calibrated D statistics support or reject a simple dosage expectation, but do not uniquely identify a backcross generation. Confirm pedigree with parental broodstock whenever possible.</div>{figures}<h2>Integrated triploid-hybrid validation</h2>{table(a.hybrid)}{cox1}<h2>Triploidy evidence</h2>{table(a.ploidy)}<h2>Global ancestry dosage</h2>{table(a.dosage)}{reciprocal}<h2>Discrete triploid dosage states and maternal consistency</h2>{table(a.states)}<h2>Window heterogeneity</h2>{table(a.windows)}<h2>Local ancestry</h2>{table(a.local)}<h2>D-statistics: excess allele sharing and calibrated P3 copies</h2>{table(a.dstat)}<h2>D-statistic synthetic mixture calibration</h2>{table(a.dcal)}{kmer}<h2>Software versions</h2>{table(a.versions)}<h2>Configuration</h2><pre>{cfg}</pre><div class='ok'>Diagnostic markers were ascertained from independent diploid parent panels, while D-statistics were calculated genome-wide and calibrated independently.</div></body></html>"""
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); tmp.write_text(body); tmp.replace(out)
if __name__=="__main__": main()
