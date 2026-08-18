#!/usr/bin/env python3
import argparse,csv,re,shutil,subprocess
from pathlib import Path
import numpy as np
import pandas as pd

def run(cmd): subprocess.run([str(x) for x in cmd],check=True)
def stat(db):
    text=subprocess.run(["meryl","statistics",str(db)],check=True,capture_output=True,text=True).stdout
    def value(key):
        m=re.search(rf"^\s*{key}\s+([0-9]+)",text,re.M|re.I); return int(m.group(1)) if m else 0
    return value("present"),value("distinct")
def recreate(path):
    if path.exists(): shutil.rmtree(path)
def input_signature(row,k,min_count):
    fields=[f"k={k}",f"min_count={min_count}"]
    for key in ("read1","read2"):
        path=Path(row[key]).resolve(); info=path.stat()
        fields.append(f"{key}={path}:{info.st_size}:{info.st_mtime_ns}")
    return "\n".join(fields)+"\n"
def pooled_fraction(g):
    a=g.parent_a_kmer_sum.sum(); return a/(a+g.parent_b_kmer_sum.sum()) if a+g.parent_b_kmer_sum.sum() else np.nan

def main():
    p=argparse.ArgumentParser(); p.add_argument("--samples",required=True); p.add_argument("--work",required=True); p.add_argument("--group-a",required=True); p.add_argument("--group-b",required=True); p.add_argument("--k",type=int,default=21); p.add_argument("--min-count",type=int,default=2); p.add_argument("--min-diagnostic-kmers",type=int,default=1000); p.add_argument("--threads",type=int,default=8); p.add_argument("--memory",type=int,default=32); p.add_argument("--output",required=True); a=p.parse_args()
    with open(a.samples) as h: rows=list(csv.DictReader(h,delimiter="\t"))
    work=Path(a.work); dbdir=work/"db"; dbdir.mkdir(parents=True,exist_ok=True); db={}
    for r in rows:
        raw=dbdir/f'{r["sample"]}.raw.meryl'; filt=dbdir/f'{r["sample"]}.filt.meryl'; done=dbdir/f'{r["sample"]}.complete'
        signature=input_signature(r,a.k,a.min_count)
        if not (filt.exists() and done.exists() and done.read_text()==signature):
            recreate(raw); recreate(filt); done.unlink(missing_ok=True)
            run(["meryl",f"k={a.k}","count",f"memory={a.memory}",f"threads={a.threads}","output",raw,r["read1"],r["read2"]]); run(["meryl","greater-than",a.min_count-1,raw,"output",filt]); shutil.rmtree(raw); done.write_text(signature)
        db[r["sample"]]=filt
    pa=[db[r["sample"]] for r in rows if r["role"]=="parent_a" and r["group"]==a.group_a]; pb=[db[r["sample"]] for r in rows if r["role"]=="parent_b" and r["group"]==a.group_b]
    if not pa or not pb: raise SystemExit("k-mer ancestry requires parent_a and parent_b panels")
    paths={name:work/f"{name}.meryl" for name in ("parent_a_core","parent_b_core","parent_a_pool","parent_b_pool","parent_a_specific","parent_b_specific")}
    for x in paths.values(): recreate(x)
    run(["meryl","intersect",*pa,"output",paths["parent_a_core"]]); run(["meryl","intersect",*pb,"output",paths["parent_b_core"]]); run(["meryl","union-sum",*pa,"output",paths["parent_a_pool"]]); run(["meryl","union-sum",*pb,"output",paths["parent_b_pool"]])
    run(["meryl","difference",paths["parent_a_core"],paths["parent_b_pool"],"output",paths["parent_a_specific"]]); run(["meryl","difference",paths["parent_b_core"],paths["parent_a_pool"],"output",paths["parent_b_specific"]])
    for label in ("parent_a_specific","parent_b_specific"):
        _,distinct=stat(paths[label])
        if distinct<a.min_diagnostic_kmers: raise SystemExit(f"{label} has only {distinct} diagnostic k-mers")
    result=[]
    for r in rows:
        vals=[]
        for label,spec in (("a",paths["parent_a_specific"]),("b",paths["parent_b_specific"])):
            tmp=work/f'tmp_{r["sample"]}_{label}.meryl'; recreate(tmp); run(["meryl","intersect",db[r["sample"]],spec,"output",tmp]); vals.append(stat(tmp)); shutil.rmtree(tmp)
        (ap,ad),(bp,bd)=vals; frac=ap/(ap+bp) if ap+bp else np.nan
        result.append({"sample":r["sample"],"group":r["group"],"role":r["role"],"ploidy":int(r["ploidy"]),"parent_a_kmer_sum":ap,"parent_b_kmer_sum":bp,"parent_a_kmer_distinct":ad,"parent_b_kmer_distinct":bd,"raw_parent_a_fraction":frac})
    frame=pd.DataFrame(result); a_signal=float(np.nanmedian(frame.loc[(frame.role=="parent_a")&(frame.group==a.group_a),"raw_parent_a_fraction"])); b_background=float(np.nanmedian(frame.loc[(frame.role=="parent_b")&(frame.group==a.group_b),"raw_parent_a_fraction"])); scale=a_signal-b_background
    if not np.isfinite(scale) or scale<=.2: raise SystemExit("k-mer parental controls do not provide a usable calibration scale")
    frame["calibrated_parent_a_fraction"]=((frame.raw_parent_a_fraction-b_background)/scale).clip(0,1); frame["calibrated_parent_a_copy_number"]=frame.calibrated_parent_a_fraction*frame.ploidy; frame["parent_a_control_signal"]=a_signal; frame["parent_b_background"]=b_background
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); frame.to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(out)

if __name__=="__main__": main()
