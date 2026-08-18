#!/usr/bin/env python3
import argparse
from pathlib import Path

def read_fasta(path):
    records=[]; name=None; seq=[]
    with open(path) as handle:
        for raw in handle:
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None: records.append((name,"".join(seq)))
                name=line[1:].split()[0]; seq=[]
            elif name is None: raise ValueError(f"sequence before FASTA header in {path}")
            else: seq.append(line.upper())
    if name is not None: records.append((name,"".join(seq)))
    return records

def main():
    p=argparse.ArgumentParser(); p.add_argument("--input",required=True); p.add_argument("--lineage",required=True); p.add_argument("--output",required=True); a=p.parse_args()
    records=read_fasta(a.input)
    if len(records)!=1: raise SystemExit(f"COX1 reference must contain exactly one sequence: {a.input}")
    seq=records[0][1]
    if not seq or set(seq)-set("ACGTN"): raise SystemExit(f"invalid COX1 sequence in {a.input}")
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp")
    with tmp.open("w") as h:
        h.write(f">{a.lineage}\n")
        for i in range(0,len(seq),80): h.write(seq[i:i+80]+"\n")
    tmp.replace(out)

if __name__=="__main__": main()
