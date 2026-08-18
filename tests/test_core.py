import importlib.util
import gzip
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
def load(name):
    p=ROOT/"workflow"/"scripts"/name; spec=importlib.util.spec_from_file_location(name,p); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def test_ancestry_summary_control_calibration_recovers_one_copy():
    m=load("summarize_ancestry.py"); rows=[]
    for sample,role,a,b,ploidy in (("T1","target",10,20,3),("A1","parent_a",30,0,2),("B1","parent_b",0,30,2)):
        for pos in range(1,5): rows.append({"sample":sample,"role":role,"tier":"strict","ploidy":ploidy,"POS":pos,"A_READS":a,"B_READS":b,"DOSAGE_DISTANCE":0.0})
    out=m.summarize(pd.DataFrame(rows)); target=out[(out["sample"]=="T1")&(out["tier"]=="strict")].iloc[0]
    assert target.calibrated_parent_a_fraction==pytest.approx(1/3); assert target.calibrated_parent_a_copy_number==pytest.approx(1.0)

def test_diagnostic_marker_requires_both_panels_near_fixed():
    m=load("diagnostic_sites.py")
    assert m.marker_orientation(.99,.01,.02)==(1,0)
    assert m.marker_orientation(.01,.99,.02)==(0,1)
    assert m.marker_orientation(.80,0.0,.10) is None

def test_dstat_components_are_reference_allele_invariant():
    m=load("dstat.py"); p=(.1,.6,.8,.05)
    assert np.allclose(m.site_components(*p),m.site_components(*(1-x for x in p)))

def test_manuscript_dstat_synthetic_calibration_recovers_copy():
    m=load("dstat.py"); p1,p3,po=.05,.95,.01
    values=[]
    for copy in np.arange(0,1.01,.1):
        p2=p1+(copy/3)*(p3-p1); x=m.site_components(p1,p2,p3,po,"manuscript_single_orientation")
        values.append({"synthetic_p3_copy_number":copy,"D":x[2]/(x[0]+x[1])})
    estimate,status=m.interpolate_copy(values[7]["D"],pd.DataFrame(values))
    assert estimate==pytest.approx(.7); assert status=="within_range"

def test_dstat_calibration_rejects_nonmonotonic_and_out_of_range():
    m=load("dstat.py")
    nonmonotonic=pd.DataFrame({"synthetic_p3_copy_number":[0,.5,1],"D":[0,.5,.4]})
    estimate,status=m.interpolate_copy(.45,nonmonotonic)
    assert np.isnan(estimate) and status=="nonmonotonic"
    monotonic=pd.DataFrame({"synthetic_p3_copy_number":[0,1,2,3],"D":[0,.2,.4,.6]})
    estimate,status=m.interpolate_copy(.7,monotonic)
    assert np.isnan(estimate) and status=="above_range"

def test_dstat_calibration_supports_decreasing_curve():
    m=load("dstat.py")
    curve=pd.DataFrame({"synthetic_p3_copy_number":[0,1,2,3],"D":[.6,.4,.2,0]})
    estimate,status=m.interpolate_copy(.3,curve)
    assert estimate==pytest.approx(1.5) and status=="within_range"

def test_dstat_null_and_positive_sign():
    m=load("dstat.py")
    assert m.site_components(.2,.2,.8,0)[2]==pytest.approx(0)
    assert m.site_components(.1,.6,.8,0)[2]>0

def test_hmm_separates_long_zero_and_two_copy_tracts():
    m=load("local_ancestry_hmm.py"); rows=[]
    for i in range(20): rows.append({"sample":"T1","CHROM":"chr1","POS":i*100000+1,"ploidy":3,"A_READS":0,"B_READS":30,"DP":30,"A_READ_FRACTION":0.0})
    for i in range(20,40): rows.append({"sample":"T1","CHROM":"chr1","POS":i*100000+1,"ploidy":3,"A_READS":20,"B_READS":10,"DP":30,"A_READ_FRACTION":2/3})
    decoded,_=m.decode(pd.DataFrame(rows),.02,[0,1,2],[0,1/3,2/3],0,1,.02)
    assert (decoded.A_COPY.iloc[:20]==0).mean()>.9; assert (decoded.A_COPY.iloc[20:]==2).mean()>.9
    segments=m.collapse(decoded,3,100000,"constrained"); assert sum(x["bp_length"] for x in segments)==decoded.POS.max()-decoded.POS.min()+1

def test_ploidy_histogram_and_lrd(tmp_path):
    m=load("summarize_ploidy.py"); hist=tmp_path/"h.tsv"; hist.write_text("33 100\n50 10\n67 100\n"); x,y=m.read_hist(hist)
    assert m.band(x,y,.31,.35)+m.band(x,y,.65,.69)==200
    lrd=tmp_path/"x.tsv"; lrd.write_text("file free dip tri tet d_dip d_tri d_tet\nx 1 1 1 1 10 2 8\n"); assert m.read_lrd(lrd)["best_model"]=="triploid"

def test_empty_ancestry_fails_cleanly():
    m=load("summarize_ancestry.py")
    with pytest.raises(ValueError): m.summarize(pd.DataFrame())

def test_input_validator_uses_supplied_references(tmp_path):
    output=tmp_path/"normalized.tsv"; ok=tmp_path/"inputs.ok"
    command=[
        sys.executable,str(ROOT/"workflow/scripts/validate_inputs.py"),
        "--reference",f"parent_a={ROOT/'tests/data/parent_a.fa'}",
        "--reference",f"parent_b={ROOT/'tests/data/parent_b.fa'}",
        "--samples",str(ROOT/"tests/data/dryrun_samples.tsv"),
        "--output",str(output),"--ok",str(ok),
    ]
    subprocess.run(command,cwd=ROOT,check=True,capture_output=True,text=True)
    assert output.is_file() and ok.read_text()=="validated\n"

def test_input_validator_accepts_sheet_without_optional_reference_column(tmp_path):
    reads1=ROOT/"tests/data/reads_1.fq"; reads2=ROOT/"tests/data/reads_2.fq"
    sheet=tmp_path/"samples.tsv"
    lines=["sample\tgroup\trole\tploidy\tread1\tread2"]
    for sample,group,role,ploidy in (("T1","T","target",3),("A1","A","parent_a",2),("B1","B","parent_b",2),("O1","O","outgroup",2)):
        lines.append(f"{sample}\t{group}\t{role}\t{ploidy}\t{reads1}\t{reads2}")
    sheet.write_text("\n".join(lines)+"\n")
    output=tmp_path/"normalized.tsv"; ok=tmp_path/"inputs.ok"
    subprocess.run([
        sys.executable,str(ROOT/"workflow/scripts/validate_inputs.py"),
        "--reference",f"parent_a={ROOT/'tests/data/parent_a.fa'}",
        "--reference",f"parent_b={ROOT/'tests/data/parent_b.fa'}",
        "--samples",str(sheet),"--output",str(output),"--ok",str(ok),
    ],cwd=ROOT,check=True,capture_output=True,text=True)
    assert "reference" not in output.read_text().splitlines()[0].split("\t")

def test_analysis_config_rejects_reused_reciprocal_reference():
    m=load("../lib.py")
    cfg={"references":{"a":"a.fa","b":"b.fa"},"analysis":{
        "primary_reference":"a","reciprocal_reference":"a",
        "parent_a_group":"A","parent_b_group":"B","target_group":"T",
        "diagnostic":{"strict_max_minor_af":.02,"moderate_max_minor_af":.1},
        "local_ancestry":{},
    }}
    with pytest.raises(ValueError,match="must differ"): m.validate_analysis_config(cfg)

def test_reference_specific_contig_selection():
    m=load("../lib.py")
    configured={"primary":["CM046776.1"],"reciprocal":["NC_088853.1"]}
    assert m.requested_contigs(configured,"primary","angulata")==["CM046776.1"]
    assert m.requested_contigs(configured,"reciprocal","gigas")==["NC_088853.1"]

def test_hybrid_validation_rejects_impossible_outgroup_requirement():
    m=load("../lib.py")
    cfg={"references":{"a":"a.fa","b":"b.fa"},"analysis":{"primary_reference":"a","reciprocal_reference":"b","parent_a_group":"A","parent_b_group":"B","target_group":"T","outgroup_group":"O","diagnostic":{"strict_max_minor_af":.02,"moderate_max_minor_af":.1},"local_ancestry":{},"hybrid_validation":{"min_significant_outgroups":2}}}
    with pytest.raises(ValueError,match="exceeds"): m.validate_analysis_config(cfg)

def test_sample_ids_reject_shell_and_path_characters(tmp_path):
    m=load("../lib.py")
    sheet=tmp_path/"samples.tsv"
    sheet.write_text("sample\tgroup\trole\tploidy\tread1\tread2\treference\nA/1\tA\tparent_a\t2\ta\tb\tr\n")
    with pytest.raises(ValueError,match="letters, numbers"):
        m.load_sample_sheet(sheet)

def test_sample_groups_must_match_configured_roles():
    m=load("../lib.py")
    cfg={"analysis":{"parent_a_group":"A","parent_b_group":"B","target_group":"T","outgroup_group":"O","outgroup_groups":["O","O2"]}}
    rows=[
        {"sample":"A1","role":"parent_a","group":"wrong","ploidy":"2"},
        {"sample":"B1","role":"parent_b","group":"B","ploidy":"2"},
        {"sample":"T1","role":"target","group":"T","ploidy":"3"},
        {"sample":"O1","role":"outgroup","group":"O","ploidy":"2"},
    ]
    with pytest.raises(ValueError,match="sample-sheet/config mismatch"):
        m.validate_sample_semantics(rows,cfg)

def test_cox1_references_require_comparable_lengths(tmp_path):
    m=load("../lib.py")
    short=tmp_path/"short.fa"; long=tmp_path/"long.fa"
    short.write_text(">short\n"+"A"*100+"\n"); long.write_text(">long\n"+"A"*150+"\n")
    with pytest.raises(ValueError,match="not comparable in length"):
        m.validate_cox1_references({"references":{"short":str(short),"long":str(long)},"max_reference_length_ratio":1.1})

def test_maternal_resolution_is_symmetric_and_discordance_is_explicit():
    m=load("dosage_states.py")
    assert m.resolve_maternal_group("A","B","A","B")== (None,"pedigree_cox1_discordant")
    assert m.resolve_maternal_group("A","A","A","B")== ("A",None)
    assert m.resolve_maternal_group(None,"B","A","B")== ("B",None)
    assert m.maternal_consistency(None,0,0,.01)=="consistent"

def test_parent_a_maternal_dosage_uses_zero_copy_as_forbidden(tmp_path):
    calls=tmp_path/"calls.tsv.gz"; cox1=tmp_path/"cox1.tsv"; output=tmp_path/"states.tsv"
    pd.DataFrame([
        {"sample":"T1","role":"target","ploidy":3,"DP":30,"CHROM":"chr1","POS":1,"NEAREST_A_COPY":0,"DOSAGE_DISTANCE":0},
        {"sample":"T1","role":"target","ploidy":3,"DP":30,"CHROM":"chr1","POS":2,"NEAREST_A_COPY":1,"DOSAGE_DISTANCE":0},
    ]).to_csv(calls,sep="\t",index=False,compression="gzip")
    pd.DataFrame([{"sample":"T1","cox1_status":"assigned","cox1_lineage":"A"}]).to_csv(cox1,sep="\t",index=False)
    subprocess.run([
        sys.executable,str(ROOT/"workflow/scripts/dosage_states.py"),"--calls",str(calls),
        "--parent-a-group","A","--parent-b-group","B","--maternal-parent-group","A",
        "--cox1-summary",str(cox1),"--output",str(output),
    ],check=True)
    row=pd.read_csv(output,sep="\t").iloc[0]
    assert row.maternal_group_evaluated=="A"
    assert row.forbidden_opposite_parent_state==0
    assert row.forbidden_opposite_parent_fraction==pytest.approx(.5)

def complete_config(**analysis_updates):
    analysis={
        "primary_reference":"r1","reciprocal_reference":"r2",
        "parent_a_group":"A","parent_b_group":"B","target_group":"T",
        "outgroup_group":"O","diagnostic":{"strict_max_minor_af":.02,"moderate_max_minor_af":.1},
        "local_ancestry":{},"dstat":{"block_size":100},"kmer":{"enabled":False},
    }
    analysis.update(analysis_updates)
    return {"references":{"r1":"a.fa","r2":"b.fa"},"analysis":analysis}

def test_triploid_target_is_required():
    m=load("../lib.py"); cfg=complete_config()
    rows=[
        {"sample":"A1","role":"parent_a","group":"A","ploidy":"2"},
        {"sample":"B1","role":"parent_b","group":"B","ploidy":"2"},
        {"sample":"T1","role":"target","group":"T","ploidy":"2"},
        {"sample":"O1","role":"outgroup","group":"O","ploidy":"2"},
    ]
    with pytest.raises(ValueError,match="ploidy=3"):
        m.validate_sample_semantics(rows,cfg)

def test_triploid_state_space_and_enabled_method_count_are_validated():
    m=load("../lib.py")
    with pytest.raises(ValueError,match="states must be within"):
        m.validate_analysis_config(complete_config(local_ancestry={"states":[0,4],"expected_fractions":[0,1]}))
    with pytest.raises(ValueError,match="exceeds enabled methods"):
        m.validate_analysis_config(complete_config(hybrid_validation={"min_methods":4}))

def test_cox1_maternal_key_and_outgroup_alias_are_validated():
    m=load("../lib.py")
    with pytest.raises(ValueError,match="COX1 reference key"):
        m.validate_analysis_config(complete_config(hybrid_validation={"maternal_parent_group":"B"},mitochondrial={"cox1":{"enabled":True,"references":{"X":"x.fa","Y":"y.fa"}}}))
    with pytest.raises(ValueError,match="included in outgroup_groups"):
        m.validate_analysis_config(complete_config(outgroup_groups=["O2"]))

def test_smoke_generator_is_small_and_has_paired_reads(tmp_path):
    output=tmp_path/"smoke"
    subprocess.run([sys.executable,str(ROOT/"tests/smoke/generate_smoke_data.py"),"--output",str(output)],check=True)
    fastqs=sorted(output.glob("*.fastq.gz"))
    assert len(fastqs)==12
    assert sum(path.stat().st_size for path in output.iterdir()) < 5_000_000
    for sample in ("T3N01","PA01","PA02","PB01","PB02","O01"):
        with gzip.open(output/f"{sample}.R1.fastq.gz","rt") as r1, gzip.open(output/f"{sample}.R2.fastq.gz","rt") as r2:
            h1=r1.readline().strip().split()[0].removesuffix("/1")
            h2=r2.readline().strip().split()[0].removesuffix("/2")
            assert h1==h2

def test_cox1_assignment_requires_identity_coverage_and_margin():
    m=load("summarize_cox1.py")
    candidates=pd.DataFrame([
        {"candidate_lineage":"ANG","score":.985,"identity_to_reference":.995,"breadth":.99},
        {"candidate_lineage":"GIG","score":.970,"identity_to_reference":.980,"breadth":.99},
    ])
    best,_,_,status=m.assign_candidate(candidates,.8,.95,.001)
    assert best.candidate_lineage=="ANG" and status=="assigned"
    _,_,_,status=m.assign_candidate(candidates,.8,.95,.02)
    assert status=="ambiguous"
