#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def fasta_sequence(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">"))

def assign_candidate(group,min_breadth,min_identity,min_margin):
    ranked=group.sort_values(["score","identity_to_reference","breadth"],ascending=False)
    best=ranked.iloc[0]; second=float(ranked.score.iloc[1]) if len(ranked)>1 else 0.0; margin=float(best.score)-second
    if best.breadth<min_breadth or best.identity_to_reference<min_identity: status="insufficient_coverage"
    elif len(ranked)>1 and margin<min_margin: status="ambiguous"
    else: status="assigned"
    return best,second,margin,status

def main():
    p=argparse.ArgumentParser(); p.add_argument("--metrics",nargs="+",required=True); p.add_argument("--min-breadth",type=float,default=.80); p.add_argument("--min-identity",type=float,default=.95); p.add_argument("--min-margin",type=float,default=.001); p.add_argument("--summary",required=True); p.add_argument("--consensus",required=True); a=p.parse_args()
    d=pd.concat([pd.read_csv(path,sep="\t") for path in a.metrics],ignore_index=True); rows=[]; selected=[]
    for sample,g in d.groupby("sample",sort=False):
        best,second,margin,status=assign_candidate(g,a.min_breadth,a.min_identity,a.min_margin)
        lineage=str(best.candidate_lineage) if status=="assigned" else status
        rows.append({"sample":sample,"cox1_lineage":lineage,"cox1_status":status,"best_candidate":best.candidate_lineage,"best_score":best.score,"second_best_score":second,"score_margin":margin,"breadth":best.breadth,"identity_to_reference":best.identity_to_reference,"mean_depth":best.mean_depth,"mapped_reads":int(best.mapped_reads),"consensus_fasta":best.consensus_fasta})
        selected.append((sample,status,str(best.candidate_lineage),fasta_sequence(best.consensus_fasta)))
    summary=Path(a.summary); summary.parent.mkdir(parents=True,exist_ok=True); tmp=summary.with_suffix(summary.suffix+".tmp"); pd.DataFrame(rows).to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(summary)
    consensus=Path(a.consensus); ctmp=consensus.with_suffix(consensus.suffix+".tmp")
    with ctmp.open("w") as h:
        for sample,status,lineage,seq in selected:
            h.write(f">{sample}|status={status}|best_candidate={lineage}\n")
            for i in range(0,len(seq),80): h.write(seq[i:i+80]+"\n")
    ctmp.replace(consensus)

if __name__=="__main__": main()
