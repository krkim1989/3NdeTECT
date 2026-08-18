#!/usr/bin/env python3
import argparse,csv,gzip
from pathlib import Path

def main():
    import pysam
    p=argparse.ArgumentParser(); p.add_argument("--vcf",required=True); p.add_argument("--sites",required=True); p.add_argument("--samples",required=True); p.add_argument("--target-group",required=True); p.add_argument("--group-a",required=True); p.add_argument("--group-b",required=True); p.add_argument("--min-dp",type=int,default=10); p.add_argument("--min-gq",type=float,default=20); p.add_argument("--output",required=True); a=p.parse_args()
    with open(a.samples) as h: meta=list(csv.DictReader(h,delimiter="\t"))
    selected={r["sample"]:r for r in meta if (r["role"]=="target" and r["group"]==a.target_group) or (r["role"]=="parent_a" and r["group"]==a.group_a) or (r["role"]=="parent_b" and r["group"]==a.group_b)}
    targets=[r for r in selected.values() if r["group"]==a.target_group and r["role"]=="target"]
    if not targets: raise SystemExit(f"no target samples for group={a.target_group}")
    sites={}
    with gzip.open(a.sites,"rt") as h:
        for r in csv.DictReader(h,delimiter="\t"): sites[(r["CHROM"],int(r["POS"]),r["REF"],r["ALT"])]=r
    if not sites: raise SystemExit("diagnostic marker file is empty")
    v=pysam.VariantFile(a.vcf); absent=sorted(set(selected)-set(v.header.samples))
    if absent: raise SystemExit("selected samples absent from VCF: "+",".join(absent))
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp")
    cols=["sample","group","role","ploidy","CHROM","POS","tier","DP","GQ","GT","A_READS","B_READS","A_READ_FRACTION","NEAREST_A_COPY","DOSAGE_DISTANCE"]
    n=0
    with gzip.open(tmp,"wt",newline="") as h:
        w=csv.DictWriter(h,fieldnames=cols,delimiter="\t",lineterminator="\n"); w.writeheader()
        for rec in v:
            if len(rec.alts or [])!=1: continue
            marker=sites.get((rec.contig,rec.pos,rec.ref,rec.alts[0]))
            if marker is None: continue
            ai=int(marker["A_ALLELE"])
            for sample,row in selected.items():
                call=rec.samples[sample]; ad=call.get("AD"); gq=call.get("GQ")
                if not ad or len(ad)!=2 or any(x is None for x in ad) or gq is None: continue
                ar=int(ad[ai]); br=int(ad[1-ai]); denom=ar+br
                if denom<a.min_dp or gq<a.min_gq: continue
                ploidy=int(row["ploidy"]); frac=ar/denom; state=max(0,min(ploidy,int(round(frac*ploidy))))
                gt="/".join("." if x is None else str(x) for x in (call.get("GT") or []))
                w.writerow({"sample":sample,"group":row["group"],"role":row["role"],"ploidy":ploidy,"CHROM":rec.contig,"POS":rec.pos,"tier":marker["tier"],"DP":denom,"GQ":gq,"GT":gt,"A_READS":ar,"B_READS":br,"A_READ_FRACTION":f"{frac:.8f}","NEAREST_A_COPY":state,"DOSAGE_DISTANCE":f"{abs(frac-state/ploidy):.8f}"}); n+=1
    if n==0: tmp.unlink(missing_ok=True); raise SystemExit("no callable ancestry observations after DP/GQ filtering")
    tmp.replace(out)

if __name__=="__main__": main()
