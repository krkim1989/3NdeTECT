#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def pooled_fraction(g):
    a=g.A_READS.sum(); total=a+g.B_READS.sum(); return a/total if total else np.nan

def summarize(df):
    required={"sample","role","tier","ploidy","POS","A_READS","B_READS","DOSAGE_DISTANCE"}
    missing=required-set(df.columns)
    if missing: raise ValueError("missing ancestry columns: "+",".join(sorted(missing)))
    rows=[]
    for tier in ("strict","moderate"):
        d=df[df.tier=="strict"] if tier=="strict" else df
        if d.empty: continue
        a_controls=d[d.role=="parent_a"]; b_controls=d[d.role=="parent_b"]
        if a_controls.empty or b_controls.empty: raise ValueError("parent controls required for ancestry calibration")
        a_signal=float(np.median([pooled_fraction(g) for _,g in a_controls.groupby("sample")]))
        b_background=float(np.median([pooled_fraction(g) for _,g in b_controls.groupby("sample")]))
        scale=a_signal-b_background
        if not np.isfinite(scale) or scale<=0.2: raise ValueError(f"invalid parental calibration scale for {tier}: {scale}")
        for sample,g in d[d.role=="target"].groupby("sample"):
            ploidy=int(g.ploidy.iloc[0]); raw=pooled_fraction(g); corrected=float(np.clip((raw-b_background)/scale,0,1))
            rows.append({"sample":sample,"tier":tier,"ploidy":ploidy,"n_sites":g.POS.nunique(),"total_allelic_depth":int((g.A_READS+g.B_READS).sum()),"raw_parent_a_fraction":raw,"parent_a_control_signal":a_signal,"parent_b_background":b_background,"calibrated_parent_a_fraction":corrected,"calibrated_parent_a_copy_number":corrected*ploidy,"pct_integer_dosage":100*(g.DOSAGE_DISTANCE<=0.10).mean()})
    if not rows: raise ValueError("no target ancestry rows available")
    return pd.DataFrame(rows)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--calls",required=True); p.add_argument("--output",required=True); a=p.parse_args()
    df=pd.read_csv(a.calls,sep="\t",compression="gzip"); result=summarize(df)
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); result.to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(out)

if __name__=="__main__": main()
