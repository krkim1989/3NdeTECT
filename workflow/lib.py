import csv
import re
from pathlib import Path

REQUIRED_SAMPLE_COLUMNS = {"sample", "group", "role", "ploidy", "read1", "read2"}
SAFE_NAME = re.compile(r"[A-Za-z0-9_.-]+")


def figure_outputs(figdir, formats, stem):
    return [f"{figdir}/{stem}.{fmt}" for fmt in formats]


def build_figure_targets(figdir, formats, enabled, painting_ancestries, kmer_enabled, tables_enabled):
    if not enabled:
        return []
    targets = []
    if "parent_a" in painting_ancestries:
        targets += figure_outputs(figdir, formats, "chromosome_painting")
    if "parent_b" in painting_ancestries:
        targets += figure_outputs(figdir, formats, "chromosome_painting_parentB")
    targets += figure_outputs(figdir, formats, "window_ancestry")
    targets += figure_outputs(figdir, formats, "dstat_calibration")
    targets += figure_outputs(figdir, formats, "integer_dosage")
    targets += figure_outputs(figdir, formats, "nquire_ploidy")
    if kmer_enabled:
        targets += figure_outputs(figdir, formats, "kmer_ancestry")
    if tables_enabled:
        targets.append(f"{figdir}/tables/manifest.tsv")
    return targets

def load_sample_sheet(path):
    with open(path, newline="") as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        missing=REQUIRED_SAMPLE_COLUMNS-set(reader.fieldnames or [])
        if missing:
            raise ValueError("missing sample-sheet columns: "+", ".join(sorted(missing)))
        rows=list(reader)
    if not rows:
        raise ValueError("sample sheet is empty")
    ids=[row["sample"].strip() for row in rows]
    if any(not sample for sample in ids) or len(ids)!=len(set(ids)):
        raise ValueError("sample IDs must be non-empty and unique")
    invalid=[sample for sample in ids if not SAFE_NAME.fullmatch(sample)]
    if invalid:
        raise ValueError("sample IDs must contain only letters, numbers, dot, underscore or hyphen: "+", ".join(invalid))
    for row,sample in zip(rows,ids):
        row["sample"]=sample
    return rows

def validate_sample_semantics(rows,cfg):
    analysis=cfg["analysis"]
    expected={
        "parent_a": {analysis["parent_a_group"]},
        "parent_b": {analysis["parent_b_group"]},
        "target": {analysis["target_group"]},
        "outgroup": set(analysis.get("outgroup_groups",[analysis["outgroup_group"]])),
    }
    errors=[]
    for role,groups in expected.items():
        members=[row for row in rows if row["role"]==role]
        if not members:
            errors.append(f"no samples with role={role}")
            continue
        wrong=[f'{row["sample"]}={row["group"]}' for row in members if row["group"] not in groups]
        if wrong:
            errors.append(f'{role} sample group must be one of {sorted(groups)}: {", ".join(wrong)}')
        represented={row["group"] for row in members}
        missing=sorted(groups-represented)
        if missing:
            errors.append(f'no {role} samples for configured group(s): {", ".join(missing)}')
    nontriploid=[]
    for row in rows:
        if row["role"]=="target":
            try:
                if int(row["ploidy"])!=3: nontriploid.append(f'{row["sample"]}={row["ploidy"]}')
            except (TypeError,ValueError):
                nontriploid.append(f'{row["sample"]}={row["ploidy"]}')
    if nontriploid:
        errors.append("3NdeTECT requires ploidy=3 for every target: "+", ".join(nontriploid))
    if errors:
        raise ValueError("sample-sheet/config mismatch:\n- "+"\n- ".join(errors))

def fasta_record_lengths(path):
    lengths=[]; current=None
    with open(path) as handle:
        for raw in handle:
            line=raw.strip()
            if not line: continue
            if line.startswith(">"):
                if current is not None: lengths.append(current)
                current=0
            elif current is None:
                raise ValueError(f"sequence before FASTA header: {path}")
            else:
                current+=len(line)
    if current is not None: lengths.append(current)
    return lengths

def validate_cox1_references(cox1):
    refs=cox1.get("references",{})
    lengths={}
    for lineage,path in refs.items():
        try:
            records=fasta_record_lengths(path)
        except OSError as exc:
            raise ValueError(f"cannot read COX1 reference {lineage}: {path}") from exc
        if len(records)!=1 or records[0]<=0:
            raise ValueError(f"COX1 reference must contain exactly one non-empty sequence: {lineage}={path}")
        lengths[lineage]=records[0]
    ratio=max(lengths.values())/min(lengths.values())
    maximum=float(cox1.get("max_reference_length_ratio",1.10))
    if ratio>maximum:
        detail=", ".join(f"{name}={length}" for name,length in lengths.items())
        raise ValueError(f"COX1 references are not comparable in length (ratio {ratio:.3f} > {maximum:.3f}): {detail}")
    return lengths

def validate_analysis_config(cfg):
    analysis=cfg["analysis"]; references=cfg["references"]
    primary=analysis["primary_reference"]; reciprocal=analysis.get("reciprocal_reference")
    if primary not in references:
        raise ValueError(f"primary_reference is not defined under references: {primary}")
    if reciprocal and reciprocal not in references:
        raise ValueError(f"reciprocal_reference is not defined under references: {reciprocal}")
    if reciprocal and reciprocal==primary:
        raise ValueError("reciprocal_reference must differ from primary_reference")
    groups=[analysis[x] for x in ("parent_a_group","parent_b_group","target_group")]
    if len(groups)!=len(set(groups)):
        raise ValueError("parent_a_group, parent_b_group and target_group must differ")
    outgroups=analysis.get("outgroup_groups",[analysis["outgroup_group"]])
    if analysis["outgroup_group"] not in outgroups:
        raise ValueError("outgroup_group must be included in outgroup_groups")
    if set(groups)&set(outgroups):
        raise ValueError("outgroup groups must differ from parent and target groups")
    hybrid=analysis.get("hybrid_validation",{})
    maternal=hybrid.get("maternal_parent_group")
    if maternal is not None and maternal not in groups[:2]:
        raise ValueError("maternal_parent_group must match parent_a_group, parent_b_group, or be null")
    if hybrid.get("min_significant_outgroups",1)>len(outgroups):
        raise ValueError("min_significant_outgroups exceeds the configured outgroup count")
    enabled_methods=3+(1 if analysis.get("kmer",{}).get("enabled",False) else 0)
    if hybrid.get("min_methods",3)>enabled_methods:
        raise ValueError(f"hybrid_validation.min_methods exceeds enabled methods ({enabled_methods})")
    expected_copy=float(hybrid.get("expected_parent_a_copies",2.0))
    if not 0<=expected_copy<=3:
        raise ValueError("expected_parent_a_copies must be within the triploid range 0..3")
    cox1=analysis.get("mitochondrial",{}).get("cox1",{})
    if cox1.get("enabled",False):
        refs=cox1.get("references",{})
        if len(refs)<2:
            raise ValueError("enabled COX1 assignment requires at least two lineage references")
        invalid=[name for name in refs if not re.fullmatch(r"[A-Za-z0-9_.-]+",name)]
        if invalid:
            raise ValueError("COX1 lineage names must contain only letters, numbers, dot, underscore or hyphen")
        if maternal is not None and maternal not in refs:
            raise ValueError("maternal_parent_group must be a COX1 reference key when COX1 is enabled")
    diagnostic=analysis["diagnostic"]
    if diagnostic["strict_max_minor_af"]>diagnostic["moderate_max_minor_af"]:
        raise ValueError("strict_max_minor_af must not exceed moderate_max_minor_af")
    local=analysis["local_ancestry"]
    states=local.get("states",[0,1,2,3]); expected=local.get("expected_fractions",[0,1/3,2/3,1])
    if len(states)!=len(expected):
        raise ValueError("local_ancestry states and expected_fractions must have equal lengths")
    if any(state<0 or state>3 for state in states):
        raise ValueError("triploid local_ancestry states must be within 0..3")
    mismatched=[(state,fraction) for state,fraction in zip(states,expected) if not abs(float(fraction)-state/3)<1e-5]
    if mismatched:
        raise ValueError("local_ancestry expected_fractions must equal state/3 for triploid targets")

def requested_contigs(configured,panel,reference_key):
    if not configured:
        return []
    if isinstance(configured,list):
        return [str(x) for x in configured]
    if not isinstance(configured,dict):
        raise ValueError("chromosomes must be a list or a panel/reference-key mapping")
    values=configured.get(panel,configured.get(reference_key,[]))
    if not isinstance(values,list):
        raise ValueError(f"chromosomes entry for {panel}/{reference_key} must be a list")
    return [str(x) for x in values]

def discover_contigs(reference, requested):
    if requested:
        return [str(x) for x in requested]
    fai=Path(f"{reference}.fai")
    if not fai.exists(): raise ValueError(f"missing FASTA index needed to discover contigs: {fai}")
    contigs=[line.split("\t",1)[0] for line in fai.read_text().splitlines() if line.strip()]
    if not contigs: raise ValueError(f"FASTA index contains no contigs: {fai}")
    return contigs

def validate_contigs(reference,requested):
    available=set(discover_contigs(reference,[]))
    absent=[contig for contig in requested if contig not in available]
    if absent:
        raise ValueError(f"contigs absent from {reference}: {','.join(absent)}")
    return requested

def reference_assets(reference):
    stem=Path(reference).with_suffix("")
    bwa=[f"{reference}.{suffix}" for suffix in ("amb","ann","bwt","pac","sa")]
    return [reference,*bwa,f"{reference}.fai",f"{stem}.dict"]
