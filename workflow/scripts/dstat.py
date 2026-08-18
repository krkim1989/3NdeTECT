#!/usr/bin/env python3
import argparse,csv,math
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

def allele_frequency(rec,names,min_called,min_dp=0,min_gq=0):
    alt=tot=0
    for sample in names:
        call=rec.samples[sample]; dp=call.get("DP"); gq=call.get("GQ")
        if dp is None or gq is None or dp<min_dp or gq<min_gq: continue
        for allele in (call.get("GT") or []):
            if allele in (0,1): tot+=1; alt+=allele
    return (alt/tot if tot>=min_called else None),tot

def site_components(p1,p2,p3,po,estimator="frequency_symmetric"):
    if estimator=="manuscript_single_orientation":
        abba=(1-p1)*p2*p3*(1-po)
        baba=p1*(1-p2)*p3*(1-po)
        numerator=abba-baba
        pd=max(p2,p3)
        fd_denominator=(1-p1)*pd*p3*(1-po)-p1*(1-pd)*p3*(1-po)
    elif estimator=="frequency_symmetric":
        abba=(1-po)*(1-p1)*p2*p3 + po*p1*(1-p2)*(1-p3)
        baba=(1-po)*p1*(1-p2)*p3 + po*(1-p1)*p2*(1-p3)
        numerator=(p2-p1)*(p3-po)
        pd_alt=max(p2,p3); pd_ref=max(1-p2,1-p3)
        fd_denominator=(1-po)*(pd_alt-p1)*pd_alt + po*(pd_ref-(1-p1))*pd_ref
    else:
        raise ValueError(f"unknown D-statistic estimator: {estimator}")
    return np.array([abba,baba,numerator,fd_denominator],dtype=float)

def ratio_and_jackknife(matrix,num_index,denominator):
    total=matrix.sum(axis=0); den=denominator(total)
    estimate=total[num_index]/den if den!=0 else np.nan; loo=[]
    for row in matrix:
        x=total-row; d=denominator(x)
        if d!=0: loo.append(x[num_index]/d)
    if len(loo)<3: return estimate,np.nan,np.nan
    arr=np.asarray(loo); se=math.sqrt((len(arr)-1)/len(arr)*np.sum((arr-arr.mean())**2)); return estimate,se,(estimate/se if se>0 else np.nan)

def bh_qvalues(pvalues):
    p=np.asarray(pvalues,float); order=np.argsort(p); q=np.full(len(p),np.nan); running=1.0
    for rank_index in range(len(p)-1,-1,-1):
        idx=order[rank_index]; running=min(running,p[idx]*len(p)/(rank_index+1)); q[idx]=running
    return q

def valid_variant(rec,min_qual):
    if len(rec.ref)!=1 or len(rec.alts or [])!=1 or len(rec.alts[0])!=1: return False
    if rec.qual is None or rec.qual<min_qual: return False
    filters=set(rec.filter.keys()); return not filters or filters=={"PASS"}

def calibration_grid(max_copy,step):
    n=int(round(max_copy/step))
    if step<=0 or max_copy<=0 or not math.isclose(n*step,max_copy,rel_tol=0,abs_tol=1e-8):
        raise ValueError("calibration max-copy must be a positive multiple of step")
    return [round(i*step,10) for i in range(n+1)]

def interpolate_copy(observed,curve):
    curve=curve[np.isfinite(curve.D)&np.isfinite(curve.synthetic_p3_copy_number)].sort_values("synthetic_p3_copy_number")
    if len(curve)<2 or not np.isfinite(observed): return np.nan,"unavailable"
    ds=curve.D.to_numpy(float); copies=curve.synthetic_p3_copy_number.to_numpy(float)
    delta=np.diff(ds); tolerance=max(1e-12,1e-9*np.nanmax(np.abs(ds)))
    increasing=np.all(delta>tolerance); decreasing=np.all(delta<-tolerance)
    if not (increasing or decreasing): return np.nan,"nonmonotonic"
    if decreasing: ds=ds[::-1]; copies=copies[::-1]
    if observed<ds[0]: return np.nan,"below_range"
    if observed>ds[-1]: return np.nan,"above_range"
    return float(np.interp(observed,ds,copies)),"within_range"

def main():
    import pysam
    p=argparse.ArgumentParser(); p.add_argument("--vcf",required=True); p.add_argument("--samples",required=True); p.add_argument("--p1",required=True); p.add_argument("--p2",required=True); p.add_argument("--p3",required=True); p.add_argument("--outgroup",action="append",required=True); p.add_argument("--block-size",type=int,default=5_000_000); p.add_argument("--min-called-alleles",type=int,default=4); p.add_argument("--min-sample-dp",type=int,default=6); p.add_argument("--min-sample-gq",type=float,default=20); p.add_argument("--min-qual",type=float,default=30); p.add_argument("--per-individual",action="store_true"); p.add_argument("--estimator",choices=["manuscript_single_orientation","frequency_symmetric"],default="manuscript_single_orientation"); p.add_argument("--calibration-max-copy",type=float,default=3.0); p.add_argument("--calibration-step",type=float,default=.1); p.add_argument("--calibration-output",required=True); p.add_argument("--output",required=True); a=p.parse_args()
    with open(a.samples) as h: meta=list(csv.DictReader(h,delimiter="\t"))
    role_by_group={a.p1:"parent_b",a.p2:"target",a.p3:"parent_a",**{group:"outgroup" for group in a.outgroup}}
    by_group={g:[r["sample"] for r in meta if r["group"]==g and r["role"]==role_by_group[g]] for g in role_by_group}
    if any(not x for x in by_group.values()): raise SystemExit("one or more D-stat groups have no samples")
    tests=[(f"{a.p2}_pooled",by_group[a.p2])]
    if a.per_individual: tests.extend((sample,[sample]) for sample in by_group[a.p2])
    grid=calibration_grid(a.calibration_max_copy,a.calibration_step)
    accum={(label,out):defaultdict(lambda:np.zeros(4)) for label,_ in tests for out in a.outgroup}; counts=defaultdict(int)
    cal={(label,out,copy):defaultdict(lambda:np.zeros(4)) for label,_ in tests for out in a.outgroup for copy in grid}
    v=pysam.VariantFile(a.vcf); required={s for members in by_group.values() for s in members}; absent=sorted(required-set(v.header.samples))
    if absent: raise SystemExit("D-stat samples absent from VCF: "+",".join(absent))
    for rec in v:
        if not valid_variant(rec,a.min_qual): continue
        p1,_=allele_frequency(rec,by_group[a.p1],a.min_called_alleles,a.min_sample_dp,a.min_sample_gq); p3,_=allele_frequency(rec,by_group[a.p3],a.min_called_alleles,a.min_sample_dp,a.min_sample_gq)
        if p1 is None or p3 is None: continue
        for label,p2names in tests:
            required_p2=1 if len(p2names)==1 else a.min_called_alleles
            p2,_=allele_frequency(rec,p2names,required_p2,a.min_sample_dp,a.min_sample_gq)
            if p2 is None: continue
            for out in a.outgroup:
                po,_=allele_frequency(rec,by_group[out],a.min_called_alleles,a.min_sample_dp,a.min_sample_gq)
                if po is None: continue
                block=(rec.contig,(rec.pos-1)//a.block_size)
                accum[(label,out)][block]+=site_components(p1,p2,p3,po,a.estimator); counts[(label,out)]+=1
                for copy in grid:
                    synthetic_p2=p1+(copy/3.0)*(p3-p1)
                    cal[(label,out,copy)][block]+=site_components(p1,synthetic_p2,p3,po,a.estimator)
    rows=[]; calrows=[]
    for (label,out),blocks in accum.items():
        if len(blocks)<3: raise SystemExit(f"{label}/{out}: fewer than three informative blocks")
        m=np.asarray(list(blocks.values())); total=m.sum(axis=0)
        D,seD,zD=ratio_and_jackknife(m,2,lambda x:x[0]+x[1]); fd,sefd,zfd=ratio_and_jackknife(m,2,lambda x:x[3])
        pval=math.erfc(abs(zD)/math.sqrt(2)) if np.isfinite(zD) else np.nan
        rows.append({"test":label,"P1":a.p1,"P2":a.p2,"P3":a.p3,"O":out,"estimator":a.estimator,"n_sites":counts[(label,out)],"n_blocks":len(m),"ABBA":total[0],"BABA":total[1],"D":D,"SE_D":seD,"Z_D":zD,"P_D":pval,"f_d":fd,"SE_fd":sefd,"Z_fd":zfd})
        for copy in grid:
            cm=np.asarray(list(cal[(label,out,copy)].values()))
            cD,cse,cz=ratio_and_jackknife(cm,2,lambda x:x[0]+x[1])
            calrows.append({"test":label,"O":out,"estimator":a.estimator,"synthetic_p3_copy_number":copy,"n_sites":counts[(label,out)],"n_blocks":len(cm),"D":cD,"SE_D":cse,"Z_D":cz})
    result=pd.DataFrame(rows); calibration=pd.DataFrame(calrows)
    estimates=[]; statuses=[]
    for row in result.itertuples(index=False):
        curve=calibration[(calibration.test==row.test)&(calibration.O==row.O)]
        estimate,status=interpolate_copy(row.D,curve); estimates.append(estimate); statuses.append(status)
    result["calibrated_p3_copy_number"]=estimates; result["calibration_status"]=statuses
    finite=result.P_D.notna(); result["Q_D_FDR"]=np.nan
    if finite.any(): result.loc[finite,"Q_D_FDR"]=bh_qvalues(result.loc[finite,"P_D"].to_numpy())
    for frame,path in ((result,Path(a.output)),(calibration,Path(a.calibration_output))):
        path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); frame.to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(path)

if __name__=="__main__": main()
