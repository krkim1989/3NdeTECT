#!/usr/bin/env python3
import argparse, csv, gzip
import re
from pathlib import Path

REQUIRED = {"sample", "group", "role", "ploidy", "read1", "read2"}
ROLES = {"target", "parent_a", "parent_b", "outgroup", "control"}
SAFE_NAME = re.compile(r"[A-Za-z0-9_.-]+")

def parse_references(items):
    refs={}
    for item in items:
        if "=" not in item: raise SystemExit(f"invalid --reference value: {item}")
        key,path=item.split("=",1)
        if not key or key in refs: raise SystemExit(f"duplicate/empty reference key: {key}")
        refs[key]=path
    return refs

def fastq_header_ok(path):
    p=Path(path)
    try:
        opener=gzip.open if p.suffix == ".gz" else open
        with opener(p,"rt") as h: return h.readline().startswith("@")
    except (OSError, UnicodeDecodeError): return False

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--reference", action="append", default=[], metavar="KEY=FASTA")
    p.add_argument("--samples", required=True); p.add_argument("--output", required=True); p.add_argument("--ok", required=True)
    a=p.parse_args(); refs=parse_references(a.reference)
    with open(a.samples, newline="") as h:
        reader=csv.DictReader(h,delimiter="\t"); fields=list(reader.fieldnames or [])
        missing=REQUIRED-set(fields)
        if missing: raise SystemExit("missing sample-sheet columns: "+", ".join(sorted(missing)))
        rows=list(reader)
    if not rows: raise SystemExit("sample sheet is empty")
    errors=[]; ids=[r["sample"].strip() for r in rows]
    if any(not x for x in ids): errors.append("sample IDs must not be empty")
    if len(ids)!=len(set(ids)): errors.append("sample IDs must be unique")
    invalid=[x for x in ids if x and not SAFE_NAME.fullmatch(x)]
    if invalid: errors.append("sample IDs must contain only letters, numbers, dot, underscore or hyphen: "+", ".join(invalid))
    for r in rows:
        sample=r["sample"]
        if r["role"] not in ROLES: errors.append(f'{sample}: invalid role {r["role"]}')
        try: ploidy=int(r["ploidy"]); assert ploidy>0
        except (ValueError, AssertionError): errors.append(f'{sample}: ploidy must be a positive integer'); ploidy=None
        if r["role"]=="target" and ploidy is not None and ploidy!=3:
            errors.append(f'{sample}: 3NdeTECT target ploidy must be 3')
        if "reference" in fields and r["reference"] and r["reference"] not in refs:
            errors.append(f'{sample}: unknown informational reference key {r["reference"]}')
        for col in ("read1","read2"):
            q=Path(r[col])
            if not q.is_file() or q.stat().st_size==0: errors.append(f'{sample}: missing/empty {col}: {q}')
            elif not fastq_header_ok(q): errors.append(f'{sample}: unreadable or invalid FASTQ header: {q}')
    for key,val in refs.items():
        fasta=Path(val)
        bwa=[Path(f"{val}.{suffix}") for suffix in ("amb","ann","bwt","pac","sa")]
        required=[fasta,*bwa,Path(val+".fai"),Path(str(fasta.with_suffix(""))+".dict")]
        for q in required:
            if not q.is_file() or q.stat().st_size==0: errors.append(f'missing/empty reference asset {key}: {q}')
    for role in ("target","parent_a","parent_b","outgroup"):
        if not any(r["role"]==role for r in rows): errors.append(f'no samples with role={role}')
    for role in ("parent_a","parent_b"):
        members=[r for r in rows if r["role"]==role]
        for r in members:
            try:
                if int(r["ploidy"])!=2: errors.append(f'{r["sample"]}: {role} marker panel must be diploid')
            except ValueError: pass
    if errors: raise SystemExit("input validation failed:\n- "+"\n- ".join(errors))
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp")
    with tmp.open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
    tmp.replace(out); ok=Path(a.ok); ok.parent.mkdir(parents=True,exist_ok=True); ok.write_text("validated\n")

if __name__=="__main__": main()
