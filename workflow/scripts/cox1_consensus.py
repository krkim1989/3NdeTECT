#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
import pysam

BASES="ACGT"

def main():
    p=argparse.ArgumentParser(); p.add_argument("--bam",required=True); p.add_argument("--reference",required=True); p.add_argument("--sample",required=True); p.add_argument("--lineage",required=True); p.add_argument("--min-depth",type=int,default=3); p.add_argument("--min-base-quality",type=int,default=20); p.add_argument("--min-major-fraction",type=float,default=.60); p.add_argument("--metrics",required=True); p.add_argument("--fasta",required=True); a=p.parse_args()
    ref=pysam.FastaFile(a.reference); contigs=list(ref.references)
    if len(contigs)!=1: raise SystemExit("normalized COX1 reference must contain one contig")
    contig=contigs[0]; sequence=ref.fetch(contig).upper(); length=len(sequence); calls=["N"]*length; depths=[0]*length
    bam=pysam.AlignmentFile(a.bam,"rb")
    for column in bam.pileup(contig,0,length,truncate=True,stepper="samtools",min_base_quality=a.min_base_quality,max_depth=100000):
        counts={base:0 for base in BASES}
        for pileup in column.pileups:
            if pileup.is_del or pileup.is_refskip: continue
            base=pileup.alignment.query_sequence[pileup.query_position].upper()
            if base in counts: counts[base]+=1
        depth=sum(counts.values()); depths[column.reference_pos]=depth
        if depth>=a.min_depth:
            base,count=max(counts.items(),key=lambda item:item[1])
            if count/depth>=a.min_major_fraction: calls[column.reference_pos]=base
    mapped=bam.count(contig=contig,read_callback="all"); bam.close(); ref.close()
    called=sum(base!="N" for base in calls); matches=sum(base!="N" and base==sequence[i] for i,base in enumerate(calls))
    breadth=called/length; identity=matches/called if called else 0.0; mean_depth=sum(depths)/length
    fasta=Path(a.fasta); fasta.parent.mkdir(parents=True,exist_ok=True); ftmp=fasta.with_suffix(fasta.suffix+".tmp")
    with ftmp.open("w") as h:
        h.write(f">{a.sample}|candidate={a.lineage}|breadth={breadth:.6f}|identity={identity:.6f}\n")
        consensus="".join(calls)
        for i in range(0,length,80): h.write(consensus[i:i+80]+"\n")
    ftmp.replace(fasta)
    metrics=Path(a.metrics); metrics.parent.mkdir(parents=True,exist_ok=True); tmp=metrics.with_suffix(metrics.suffix+".tmp")
    with tmp.open("w",newline="") as h:
        fields=["sample","candidate_lineage","reference_length","mapped_reads","called_bases","breadth","identity_to_reference","mean_depth","score","consensus_fasta"]
        w=csv.DictWriter(h,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerow({"sample":a.sample,"candidate_lineage":a.lineage,"reference_length":length,"mapped_reads":mapped,"called_bases":called,"breadth":f"{breadth:.7g}","identity_to_reference":f"{identity:.7g}","mean_depth":f"{mean_depth:.7g}","score":f"{breadth*identity:.7g}","consensus_fasta":str(fasta)})
    tmp.replace(metrics)

if __name__=="__main__": main()
