#!/usr/bin/env python3
import argparse, csv, gzip
from pathlib import Path

def panel_af(record, names, min_dp=0, min_gq=0):
    alt=called=0
    for name in names:
        call=record.samples[name]
        gt=call.get("GT")
        dp=call.get("DP"); gq=call.get("GQ")
        if not gt or dp is None or gq is None or dp<min_dp or gq<min_gq: continue
        for allele in gt:
            if allele in (0,1): called+=1; alt+=allele
    return (alt/called if called else None),called

def marker_orientation(af_a,af_b,threshold):
    if af_a>=1-threshold and af_b<=threshold: return 1,0
    if af_a<=threshold and af_b>=1-threshold: return 0,1
    return None

def usable_record(rec,min_qual):
    if len(rec.ref)!=1 or len(rec.alts or [])!=1 or len(rec.alts[0])!=1: return False
    if rec.qual is None or rec.qual<min_qual: return False
    filters=set(rec.filter.keys())
    return not filters or filters=={"PASS"}

def main():
    import pysam
    p=argparse.ArgumentParser(); p.add_argument("--vcf",required=True); p.add_argument("--samples",required=True)
    p.add_argument("--group-a",required=True); p.add_argument("--group-b",required=True)
    p.add_argument("--strict",type=float,default=.02); p.add_argument("--moderate",type=float,default=.10)
    p.add_argument("--min-called-alleles",type=int,default=4); p.add_argument("--max-missing-fraction",type=float,default=.2)
    p.add_argument("--min-parent-dp",type=int,default=6); p.add_argument("--min-parent-gq",type=float,default=20); p.add_argument("--min-qual",type=float,default=30)
    p.add_argument("--output",required=True); a=p.parse_args()
    if not 0<=a.strict<=a.moderate<.5: raise SystemExit("require 0 <= strict <= moderate < 0.5")
    with open(a.samples) as h: meta=list(csv.DictReader(h,delimiter="\t"))
    ga=[r["sample"] for r in meta if r["group"]==a.group_a and r["role"]=="parent_a"]
    gb=[r["sample"] for r in meta if r["group"]==a.group_b and r["role"]=="parent_b"]
    if not ga or not gb: raise SystemExit("parent groups not found in independent diploid panels")
    v=pysam.VariantFile(a.vcf); missing_samples=sorted((set(ga)|set(gb))-set(v.header.samples))
    if missing_samples: raise SystemExit("parent samples absent from VCF: "+",".join(missing_samples))
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp")
    cols=["CHROM","POS","REF","ALT","QUAL","A_ALT_AF","B_ALT_AF","A_N","B_N","A_MISSING","B_MISSING","A_ALLELE","B_ALLELE","tier"]
    n_written=0
    with gzip.open(tmp,"wt",newline="") as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter="\t",lineterminator="\n"); w.writeheader()
        for rec in v:
            if not usable_record(rec,a.min_qual): continue
            af_a,n_a=panel_af(rec,ga,a.min_parent_dp,a.min_parent_gq); af_b,n_b=panel_af(rec,gb,a.min_parent_dp,a.min_parent_gq)
            miss_a=1-n_a/(2*len(ga)); miss_b=1-n_b/(2*len(gb))
            if af_a is None or af_b is None or min(n_a,n_b)<a.min_called_alleles or max(miss_a,miss_b)>a.max_missing_fraction: continue
            orient=marker_orientation(af_a,af_b,a.strict); tier="strict"
            if orient is None: orient=marker_orientation(af_a,af_b,a.moderate); tier="moderate"
            if orient is None: continue
            ai,bi=orient
            w.writerow({"CHROM":rec.contig,"POS":rec.pos,"REF":rec.ref,"ALT":rec.alts[0],"QUAL":f"{rec.qual:.3f}","A_ALT_AF":f"{af_a:.6f}","B_ALT_AF":f"{af_b:.6f}","A_N":n_a,"B_N":n_b,"A_MISSING":f"{miss_a:.6f}","B_MISSING":f"{miss_b:.6f}","A_ALLELE":ai,"B_ALLELE":bi,"tier":tier}); n_written+=1
    if n_written==0: tmp.unlink(missing_ok=True); raise SystemExit("no diagnostic markers passed all filters")
    tmp.replace(out)

if __name__=="__main__": main()
