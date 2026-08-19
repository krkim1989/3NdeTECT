#!/usr/bin/env python3
import argparse,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import betaln,gammaln

def log_beta_binom(k,n,mu,rho):
    mu=float(np.clip(mu,1e-5,1-1e-5)); rho=float(np.clip(rho,1e-6,.95)); concentration=(1-rho)/rho
    alpha=mu*concentration; beta=(1-mu)*concentration
    return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+betaln(k+alpha,n-k+beta)-betaln(alpha,beta)

def decode(group,transition_per_mb,states,expected,b_background=0.0,a_signal=1.0,rho=.02):
    g=group.sort_values("POS").reset_index(drop=True)
    if g.empty: raise ValueError("cannot decode an empty chromosome")
    means=b_background+(a_signal-b_background)*np.asarray(expected,dtype=float)
    n=len(g); m=len(states); score=np.full((n,m),-np.inf); back=np.zeros((n,m),dtype=int)
    for j,mu in enumerate(means): score[0,j]=-math.log(m)+log_beta_binom(int(g.A_READS.iloc[0]),int(g.DP.iloc[0]),mu,rho)
    for i in range(1,n):
        dist=max(0,int(g.POS.iloc[i])-int(g.POS.iloc[i-1]))/1e6
        switch=min(.49,1-math.exp(-transition_per_mb*dist)); stay=max(1-switch,1e-12); other=max(switch/max(m-1,1),1e-12)
        for j,mu in enumerate(means):
            z=score[i-1]+np.array([math.log(stay if k==j else other) for k in range(m)])
            back[i,j]=int(np.argmax(z)); score[i,j]=z[back[i,j]]+log_beta_binom(int(g.A_READS.iloc[i]),int(g.DP.iloc[i]),mu,rho)
    path=np.zeros(n,dtype=int); path[-1]=int(np.argmax(score[-1]))
    for i in range(n-2,-1,-1): path[i]=back[i+1,path[i+1]]
    g["A_COPY"]=[states[x] for x in path]; return g,float(np.max(score[-1]))

def collapse(g,min_sites,min_bp,model):
    rows=[]; start=0; positions=g.POS.to_numpy(dtype=int)
    for i in range(1,len(g)+1):
        if i==len(g) or g.A_COPY.iloc[i]!=g.A_COPY.iloc[start]:
            last=i-1; left=positions[0] if start==0 else (positions[start-1]+positions[start])//2+1
            right=positions[-1] if last==len(g)-1 else (positions[last]+positions[last+1])//2
            z=g.iloc[start:i]; bp=max(1,right-left+1)
            rows.append({"sample":z["sample"].iloc[0],"model":model,"CHROM":z["CHROM"].iloc[0],"start":left,"end":right,"bp_length":bp,"ploidy":int(z["ploidy"].iloc[0]),"A_COPY":int(z.A_COPY.iloc[0]),"n_sites":len(z),"mean_A_read_fraction":z.A_READ_FRACTION.mean(),"pass_filter":len(z)>=min_sites and bp>=min_bp})
            start=i
    return rows

def pooled(g):
    a=g.A_READS.sum(); return a/(a+g.B_READS.sum())

def main():
    p=argparse.ArgumentParser(); p.add_argument("--calls",required=True); p.add_argument("--transition-per-mb",type=float,default=.015); p.add_argument("--states",default="0,1,2,3"); p.add_argument("--expected-fractions",default="0,0.333333,0.666667,1"); p.add_argument("--rho",type=float,default=.02); p.add_argument("--min-sites",type=int,default=3); p.add_argument("--min-bp",type=int,default=1_000_000); p.add_argument("--compare-unconstrained",action="store_true"); p.add_argument("--segments",required=True); p.add_argument("--summary",required=True); a=p.parse_args()
    states=[int(x) for x in a.states.split(",")]; expected=[float(x) for x in a.expected_fractions.split(",")]
    if len(states)!=len(expected): raise SystemExit("states and expected-fractions lengths differ")
    df=pd.read_csv(a.calls,sep="\t",compression="gzip"); targets=df[df.role=="target"].copy()
    if targets.empty: raise SystemExit("no target calls for local ancestry")
    a_controls=df[df.role=="parent_a"]; b_controls=df[df.role=="parent_b"]
    if a_controls.empty or b_controls.empty: raise SystemExit("parent controls required for HMM emission calibration")
    a_signal=float(np.median([pooled(g) for _,g in a_controls.groupby("sample")]))
    b_background=float(np.median([pooled(g) for _,g in b_controls.groupby("sample")]))
    if a_signal-b_background<=.2: raise SystemExit("parental HMM calibration has insufficient separation")
    seg=[]; scores={}
    for (sample,chrom),g in targets.groupby(["sample","CHROM"],sort=False):
        d,score=decode(g,a.transition_per_mb,states,expected,b_background,a_signal,a.rho); seg.extend(collapse(d,a.min_sites,a.min_bp,"constrained")); scores[(sample,chrom,"constrained")]=score
        if a.compare_unconstrained:
            ploidy=int(g.ploidy.iloc[0]); full_states=list(range(ploidy+1)); full_expected=[x/ploidy for x in full_states]
            u,uscore=decode(g,a.transition_per_mb,full_states,full_expected,b_background,a_signal,a.rho); seg.extend(collapse(u,a.min_sites,a.min_bp,"unconstrained")); scores[(sample,chrom,"unconstrained")]=uscore
    ss=pd.DataFrame(seg)
    if ss.empty: raise SystemExit("HMM produced no segments")
    out=Path(a.segments); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); ss.to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(out)
    rows=[]
    for sample,g in ss[ss.model=="constrained"].groupby("sample"):
        ploidy=int(targets[targets["sample"]==sample].ploidy.iloc[0]); total=g.bp_length.sum(); passed=g[g["pass_filter"].astype(bool)]; pass_bp=passed.bp_length.sum()
        full=ss[(ss.sample==sample)&(ss.model=="unconstrained")]; full_total=full.bp_length.sum()
        rows.append({"sample":sample,"ploidy":ploidy,"n_segments":len(g),"n_pass_segments":len(passed),"marker_span_bp":int(total),"bpweighted_parent_a_ancestry":float((g.A_COPY*g.bp_length).sum()/(ploidy*total)),"filtered_bpweighted_parent_a_ancestry":float((passed.A_COPY*passed.bp_length).sum()/(ploidy*pass_bp)) if pass_bp else np.nan,"unconstrained_3copy_bp_fraction":float(full.loc[full.A_COPY==ploidy,"bp_length"].sum()/full_total) if full_total else np.nan,"viterbi_loglik_constrained":sum(v for (s,c,m),v in scores.items() if s==sample and m=="constrained"),"viterbi_loglik_unconstrained":sum(v for (s,c,m),v in scores.items() if s==sample and m=="unconstrained") if a.compare_unconstrained else np.nan,"parent_a_control_signal":a_signal,"parent_b_background":b_background,"rho":a.rho})
    summary=Path(a.summary); tmp=summary.with_suffix(summary.suffix+".tmp"); pd.DataFrame(rows).to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(summary)

if __name__=="__main__": main()
