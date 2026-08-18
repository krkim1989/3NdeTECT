#!/usr/bin/env python3
import argparse
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd

def longest_run(states,value):
    best=current=0
    for state in states:
        current=current+1 if state==value else 0; best=max(best,current)
    return best

def resolve_maternal_group(pedigree,cox_lineage,parent_a,parent_b):
    if pedigree and cox_lineage and pedigree!=cox_lineage:
        return None,"pedigree_cox1_discordant"
    resolved=pedigree or cox_lineage
    if resolved is None: return None,"not_evaluated"
    if resolved not in {parent_a,parent_b}: return None,"cox1_nonparent_lineage"
    return resolved,None

def maternal_consistency(error_status,forbidden_fraction,max_run,max_fraction):
    if error_status: return error_status
    return "consistent" if forbidden_fraction<=max_fraction and max_run<=1 else "review"

def main():
    p=argparse.ArgumentParser(); p.add_argument("--calls",required=True); p.add_argument("--min-dp",type=int,default=20); p.add_argument("--tolerance",type=float,default=.10); p.add_argument("--max-forbidden-fraction",type=float,default=.01); p.add_argument("--maternal-parent-group"); p.add_argument("--cox1-summary"); p.add_argument("--parent-a-group",required=True); p.add_argument("--parent-b-group",required=True); p.add_argument("--output",required=True); a=p.parse_args()
    d=pd.read_csv(a.calls,sep="\t",compression="gzip"); d=d[(d.role=="target")&(d.DP>=a.min_dp)].copy()
    if d.empty: raise SystemExit("no high-depth target calls for dosage-state analysis")
    d["site_key"]=d.CHROM.astype(str)+":"+d.POS.astype(str)
    cox1=pd.read_csv(a.cox1_summary,sep="\t").set_index("sample") if a.cox1_summary else None
    prepared=[]; all_forbidden=[]
    for sample,g in d.groupby("sample",sort=False):
        ploidy=int(g.ploidy.iloc[0]); counts=g.NEAREST_A_COPY.value_counts(); n=len(g)
        row={"sample":sample,"ploidy":ploidy,"min_dp":a.min_dp,"n_sites":n,"pct_near_integer_dosage":100*(g.DOSAGE_DISTANCE<=a.tolerance).mean()}
        for state in range(ploidy+1): row[f"state_{state}_sites"]=int(counts.get(state,0)); row[f"state_{state}_fraction"]=counts.get(state,0)/n
        cox_assigned=cox1 is not None and sample in cox1.index and cox1.loc[sample,"cox1_status"]=="assigned"
        cox_lineage=str(cox1.loc[sample,"cox1_lineage"]) if cox_assigned else None
        resolved,error=resolve_maternal_group(a.maternal_parent_group,cox_lineage,a.parent_a_group,a.parent_b_group)
        forbidden_state=0 if resolved==a.parent_a_group else ploidy if resolved==a.parent_b_group else None
        if forbidden_state is None:
            forbidden_sites=g.iloc[0:0]; ff=np.nan; maxrun=np.nan
        else:
            forbidden_sites=g[g.NEAREST_A_COPY==forbidden_state]; ff=len(forbidden_sites)/n
            ordered=g.sort_values(["CHROM","POS"])
            maxrun=max((longest_run(c.NEAREST_A_COPY.tolist(),forbidden_state) for _,c in ordered.groupby("CHROM")),default=0)
        maternal=maternal_consistency(error,ff,maxrun,a.max_forbidden_fraction) if forbidden_state is not None else error
        keys=[(key,forbidden_state) for key in forbidden_sites.site_key]; all_forbidden.extend(keys)
        row.update({"maternal_group_pedigree":a.maternal_parent_group or "not_provided","maternal_group_cox1":cox_lineage or "not_assigned","maternal_group_evaluated":resolved or "not_evaluated","forbidden_opposite_parent_state":forbidden_state if forbidden_state is not None else np.nan,"forbidden_opposite_parent_fraction":ff,"max_forbidden_run_sites":maxrun,"maternal_lineage_consistency":maternal})
        prepared.append((row,keys))
    recurrence=Counter(all_forbidden); rows=[]
    for row,keys in prepared:
        row["recurrent_forbidden_sites"]=int(sum(recurrence[key]>=2 for key in keys)); rows.append(row)
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); pd.DataFrame(rows).to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(out)

if __name__=="__main__": main()
