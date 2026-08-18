#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def main():
    p=argparse.ArgumentParser(); p.add_argument("--primary",required=True); p.add_argument("--reciprocal",required=True); p.add_argument("--output",required=True); a=p.parse_args()
    x=pd.read_csv(a.primary,sep="\t"); y=pd.read_csv(a.reciprocal,sep="\t")
    z=x.merge(y,on=["sample","tier","ploidy"],suffixes=("_primary","_reciprocal"))
    if z.empty: raise SystemExit("no matching primary/reciprocal ancestry rows")
    z["absolute_copy_difference"]=(z.calibrated_parent_a_copy_number_primary-z.calibrated_parent_a_copy_number_reciprocal).abs()
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); z.to_csv(out,sep="\t",index=False,float_format="%.6f")

if __name__=="__main__": main()
