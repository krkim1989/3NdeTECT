#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def pooled(g):
    a=g.A_READS.sum(); return a/(a+g.B_READS.sum())

def main():
    p=argparse.ArgumentParser(); p.add_argument("--calls",required=True); p.add_argument("--window-size",type=int,default=5_000_000); p.add_argument("--min-sites",type=int,default=10); p.add_argument("--windows",required=True); p.add_argument("--summary",required=True); a=p.parse_args()
    d=pd.read_csv(a.calls,sep="\t",compression="gzip")
    targets=d[d.role=="target"].copy(); ac=d[d.role=="parent_a"]; bc=d[d.role=="parent_b"]
    if targets.empty or ac.empty or bc.empty: raise SystemExit("target and both parent controls are required")
    a_signal=float(np.median([pooled(g) for _,g in ac.groupby("sample")])); b_background=float(np.median([pooled(g) for _,g in bc.groupby("sample")]))
    scale=a_signal-b_background
    if scale<=.2: raise SystemExit("insufficient parental separation for window calibration")
    targets["WIN_START"]=(targets.POS//a.window_size)*a.window_size
    w=(targets.groupby(["sample","CHROM","WIN_START"]).agg(A_READS=("A_READS","sum"),B_READS=("B_READS","sum"),n_sites=("POS","nunique"),ploidy=("ploidy","first")).reset_index())
    w=w[w.n_sites>=a.min_sites].copy()
    if w.empty: raise SystemExit("no windows passed the minimum marker count")
    w["raw_parent_a_fraction"]=w.A_READS/(w.A_READS+w.B_READS); w["calibrated_parent_a_fraction"]=((w.raw_parent_a_fraction-b_background)/scale).clip(0,1); w["calibrated_parent_a_copy_number"]=w.calibrated_parent_a_fraction*w.ploidy
    out=Path(a.windows); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); w.to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(out)
    rows=[]
    for sample,g in w.groupby("sample"):
        depths=g.A_READS+g.B_READS; pbar=np.average(g.raw_parent_a_fraction,weights=depths); obs=g.raw_parent_a_fraction.std(ddof=1) if len(g)>1 else np.nan; expected=np.sqrt(np.mean(pbar*(1-pbar)/depths))
        rows.append({"sample":sample,"n_windows":len(g),"raw_parent_a_fraction":pbar,"calibrated_parent_a_fraction":float(np.clip((pbar-b_background)/scale,0,1)),"sd_across_windows":obs,"expected_binomial_sd":expected,"overdispersion_ratio":obs/expected if expected and np.isfinite(obs) else np.nan,"window_min":g.calibrated_parent_a_fraction.min(),"window_max":g.calibrated_parent_a_fraction.max(),"parent_a_control_signal":a_signal,"parent_b_background":b_background})
    summary=Path(a.summary); tmp=summary.with_suffix(summary.suffix+".tmp"); pd.DataFrame(rows).to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(summary)

if __name__=="__main__": main()
