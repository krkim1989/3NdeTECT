#!/usr/bin/env python3
"""Integrate independent triploidy and hybrid-ancestry evidence per target."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def read(path): return pd.read_csv(path,sep="\t")
def sample_row(frame,sample,**filters):
    x=frame[frame["sample"]==sample]
    for key,value in filters.items(): x=x[x[key]==value]
    return x.iloc[0] if not x.empty else None
def finite(value): return value is not None and np.isfinite(float(value))

def main():
    p=argparse.ArgumentParser()
    for name in ("ploidy","dosage","states","local","dstat"): p.add_argument(f"--{name}",required=True)
    p.add_argument("--reciprocal"); p.add_argument("--kmer"); p.add_argument("--cox1"); p.add_argument("--tier",default="strict"); p.add_argument("--min-parent-a-copy",type=float,default=.10); p.add_argument("--min-methods",type=int,default=3); p.add_argument("--min-significant-outgroups",type=int,default=1); p.add_argument("--alpha",type=float,default=.05); p.add_argument("--expected-parent-a-copies",type=float,default=2.0); p.add_argument("--simple-cross-tolerance",type=float,default=.35); p.add_argument("--output",required=True)
    a=p.parse_args(); ploidy=read(a.ploidy); dosage=read(a.dosage); states=read(a.states); local=read(a.local); dstat=read(a.dstat); reciprocal=read(a.reciprocal) if a.reciprocal else None; kmer=read(a.kmer) if a.kmer else None; cox1=read(a.cox1) if a.cox1 else None
    targets=dosage[(dosage.tier==a.tier)].copy()
    if targets.empty: raise SystemExit(f"no target dosage rows for tier {a.tier}")
    rows=[]
    for x in targets.itertuples(index=False):
        sample=x.sample; ncopy=float(x.calibrated_parent_a_copy_number); nploidy=int(x.ploidy)
        prow=sample_row(ploidy,sample); srow=sample_row(states,sample); lrow=sample_row(local,sample)
        drows=dstat[dstat.test==sample]
        if drows.empty: drows=dstat[dstat.test.str.endswith("_pooled")]
        significant=drows[(drows.D>0)&(drows.Q_D_FDR<=a.alpha)]
        if "calibration_status" in drows and not drows.empty:
            dstatuses=sorted(set(drows.calibration_status.astype(str)))
            valid_drows=drows[drows.calibration_status=="within_range"]
            dstatus="within_range" if len(valid_drows)==len(drows) else ";".join(dstatuses)
        else:
            valid_drows=drows
            dstatus="legacy_or_unavailable"
        dcopy=float(valid_drows.calibrated_p3_copy_number.mean()) if "calibrated_p3_copy_number" in valid_drows and not valid_drows.empty else np.nan
        rrow=sample_row(reciprocal,sample,tier=a.tier) if reciprocal is not None else None
        krow=sample_row(kmer,sample) if kmer is not None else None
        crow=sample_row(cox1,sample) if cox1 is not None else None
        local_copy=float(lrow.bpweighted_parent_a_ancestry*nploidy) if lrow is not None else np.nan
        kcopy=float(krow.calibrated_parent_a_copy_number) if krow is not None else np.nan
        dosage_support=a.min_parent_a_copy<ncopy<nploidy-a.min_parent_a_copy
        d_support=len(significant)>=a.min_significant_outgroups
        local_support=finite(local_copy) and a.min_parent_a_copy<local_copy<nploidy-a.min_parent_a_copy
        kmer_support=finite(kcopy) and a.min_parent_a_copy<kcopy<nploidy-a.min_parent_a_copy
        methods={"dosage":dosage_support,"dstat":d_support,"local_ancestry":local_support}
        if kmer is not None: methods["kmer"]=kmer_support
        evidence_count=sum(methods.values()); ploidy_call=prow.ploidy_call if prow is not None else "missing"
        triploid_ok=ploidy_call in {"triploid_supported","triploid_probable"}
        hybrid_call=("triploid_hybrid_supported" if triploid_ok and evidence_count>=a.min_methods else "hybrid_supported_ploidy_unresolved" if evidence_count>=a.min_methods else "hybrid_not_resolved")
        simple=abs(ncopy-a.expected_parent_a_copies)<=a.simple_cross_tolerance
        nstates=sum(int(srow.get(f"state_{i}_sites",0))>0 for i in range(nploidy+1)) if srow is not None else 0
        cox_lineage=crow.cox1_lineage if crow is not None else "not_evaluated"; cox_status=crow.cox1_status if crow is not None else "not_evaluated"
        rows.append({"sample":sample,"ploidy_call":ploidy_call,"hybrid_call":hybrid_call,"independent_hybrid_methods_supported":evidence_count,"methods_evaluated":len(methods),"dosage_support":dosage_support,"dstat_support":d_support,"local_ancestry_support":local_support,"kmer_support":kmer_support if kmer is not None else np.nan,"parent_a_copy_dosage":ncopy,"parent_a_copy_dstat_calibrated":dcopy,"dstat_calibration_status":dstatus,"parent_a_copy_kmer":kcopy,"parent_a_copy_local_ancestry":local_copy,"significant_outgroups":len(significant),"outgroups_evaluated":drows.O.nunique(),"reciprocal_absolute_copy_difference":float(rrow.absolute_copy_difference) if rrow is not None else np.nan,"pct_integer_dosage":float(x.pct_integer_dosage),"n_observed_dosage_states":nstates,"cox1_lineage":cox_lineage,"cox1_status":cox_status,"cox1_matches_configured_maternal_group":cox_status=="assigned" and srow is not None and srow.maternal_group_pedigree!="not_provided" and cox_lineage==srow.maternal_group_pedigree,"maternal_group_evaluated":srow.maternal_group_evaluated if srow is not None else "not_evaluated","forbidden_opposite_parent_state":float(srow.forbidden_opposite_parent_state) if srow is not None and finite(srow.forbidden_opposite_parent_state) else np.nan,"forbidden_opposite_parent_fraction":float(srow.forbidden_opposite_parent_fraction) if srow is not None and finite(srow.forbidden_opposite_parent_fraction) else np.nan,"max_forbidden_run_sites":int(srow.max_forbidden_run_sites) if srow is not None and finite(srow.max_forbidden_run_sites) else np.nan,"maternal_lineage_consistency":srow.maternal_lineage_consistency if srow is not None else "not_evaluated","expected_simple_cross_parent_a_copies":a.expected_parent_a_copies,"simple_cross_compatible":simple,"pedigree_interpretation":"compatible_with_simple_cross_dosage" if simple else "inconsistent_with_simple_cross_dosage" if hybrid_call!="hybrid_not_resolved" else "insufficient_hybrid_evidence"})
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); pd.DataFrame(rows).to_csv(tmp,sep="\t",index=False,float_format="%.7g"); tmp.replace(out)

if __name__=="__main__": main()
