#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
import numpy as np
import pandas as pd

def read_hist(path):
    vals=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"): continue
        f=line.split()
        try: vals.append((float(f[0]),float(f[1])))
        except (ValueError,IndexError): continue
    if not vals: raise ValueError(f"no numeric histogram rows in {path}")
    x=np.array([z[0] for z in vals]); y=np.array([z[1] for z in vals]); x=x/100 if x.max()>1.5 else x; return x,y
def band(x,y,lo,hi): return float(y[(x>=lo)&(x<=hi)].sum())
def read_lrd(path):
    lines=[x.split() for x in Path(path).read_text().splitlines() if x.strip()]
    if len(lines)<2 or len(lines[-1])<8: raise ValueError(f"invalid lrdmodel output: {path}")
    vals=[float(lines[-1][5]),float(lines[-1][6]),float(lines[-1][7])]; names=["diploid","triploid","tetraploid"]
    return {"delta_diploid":vals[0],"delta_triploid":vals[1],"delta_tetraploid":vals[2],"best_model":names[int(np.argmin(vals))]}
def classify(row,min_density,threshold):
    if row["informative_density"]<min_density: return "insufficient_data"
    model=row["best_model"]=="triploid"; spectrum=row["thirds_to_half_ratio"]>threshold
    if model and spectrum: return "triploid_supported"
    if model or spectrum: return "triploid_probable"
    return "triploid_not_supported" if row["expected_ploidy"]==3 else "nontriploid_control"
def main():
    p=argparse.ArgumentParser(); p.add_argument("--samples",required=True); p.add_argument("--hist-dir",required=True); p.add_argument("--min-spectrum-density",type=float,default=1000); p.add_argument("--diploid-ratio-quantile",type=float,default=.95); p.add_argument("--mapping-panel",default="primary"); p.add_argument("--output",required=True); a=p.parse_args()
    with open(a.samples) as h: samples=list(csv.DictReader(h,delimiter="\t"))
    rows=[]
    for s in samples:
        x,y=read_hist(Path(a.hist_dir)/f'{s["sample"]}.hist.tsv'); thirds=band(x,y,.31,.35)+band(x,y,.65,.69); half=band(x,y,.48,.52); model=read_lrd(Path(a.hist_dir)/f'{s["sample"]}.lrd.denoised.tsv')
        rows.append({"sample":s["sample"],"group":s["group"],"role":s["role"],"mapping_panel":a.mapping_panel,"expected_ploidy":int(s["ploidy"]),**model,"thirds_density":thirds,"half_density":half,"informative_density":thirds+half,"thirds_to_half_ratio":thirds/max(half,1e-12)})
    frame=pd.DataFrame(rows); thresholds={}
    for panel,g in frame[(frame.expected_ploidy==2)&(frame.role!="target")].groupby("mapping_panel"): thresholds[panel]=float(g.thirds_to_half_ratio.quantile(a.diploid_ratio_quantile))
    if set(frame.mapping_panel)-set(thresholds): raise SystemExit("each mapping panel requires at least one diploid control")
    frame["diploid_ratio_threshold"]=frame.mapping_panel.map(thresholds); frame["ploidy_call"]=[classify(r,a.min_spectrum_density,r["diploid_ratio_threshold"]) for _,r in frame.iterrows()]
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); frame.to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(out)
if __name__=="__main__": main()
