#!/usr/bin/env python3
"""Generate a deterministic, small triploid-hybrid whole-workflow dataset."""

import argparse
import gzip
import io
import random
from pathlib import Path

LENGTH = 200_000
READ_LENGTH = 100
FRAGMENT = 260
STEP = 50
COMPLEMENT = str.maketrans("ACGT", "TGCA")
ALT = {"A": "C", "C": "G", "G": "T", "T": "A"}


def mutate(sequence, positions):
    bases=list(sequence)
    for position in positions:
        bases[position]=ALT[bases[position]]
    return "".join(bases)


def write_fasta(path, contig, sequence):
    with path.open("w") as handle:
        handle.write(f">{contig}\n")
        for start in range(0,len(sequence),80):
            handle.write(sequence[start:start+80]+"\n")


def gzip_text(path):
    raw=path.open("wb")
    stream=gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0)
    return io.TextIOWrapper(stream,encoding="ascii",newline="\n")


def write_reads(path1,path2,sample,haplotypes):
    quality="I"*READ_LENGTH
    with gzip_text(path1) as r1, gzip_text(path2) as r2:
        for copy,sequence in enumerate(haplotypes):
            offset=(copy*11+len(sample)*7)%STEP
            for start in range(offset,len(sequence)-FRAGMENT,STEP):
                name=f"{sample}_{copy}_{start}"
                forward=sequence[start:start+READ_LENGTH]
                reverse=sequence[start+FRAGMENT-READ_LENGTH:start+FRAGMENT].translate(COMPLEMENT)[::-1]
                r1.write(f"@{name}/1\n{forward}\n+\n{quality}\n")
                r2.write(f"@{name}/2\n{reverse}\n+\n{quality}\n")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="tests/smoke/data")
    args=parser.parse_args()
    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    rng=random.Random(3_141_592)
    reference="".join(rng.choice("ACGT") for _ in range(LENGTH))
    common=set(range(503,LENGTH-503,211))
    diagnostic=set(range(1009,LENGTH-503,307))-common
    common_alt=mutate(reference,common)
    parent_a_ref=mutate(reference,diagnostic)
    parent_a_alt=mutate(parent_a_ref,common)

    write_fasta(out/"primary.fa","pri_chr1",reference)
    write_fasta(out/"reciprocal.fa","rec_chr1",reference)
    cox1_start,cox1_end=50_000,51_200
    write_fasta(out/"parent_a_cox1.fa","parent_a",parent_a_ref[cox1_start:cox1_end])
    write_fasta(out/"parent_b_cox1.fa","parent_b",reference[cox1_start:cox1_end])

    samples={
        "T3N01":[parent_a_ref,reference,common_alt],
        "PA01":[parent_a_ref,parent_a_alt],
        "PA02":[parent_a_alt,parent_a_ref],
        "PB01":[reference,common_alt],
        "PB02":[common_alt,reference],
        "O01":[reference,common_alt],
    }
    for sample,haplotypes in samples.items():
        write_reads(out/f"{sample}.R1.fastq.gz",out/f"{sample}.R2.fastq.gz",sample,haplotypes)
    (out/"truth.tsv").write_text(
        "quantity\tvalue\n"
        f"reference_length\t{LENGTH}\n"
        f"diagnostic_sites\t{len(diagnostic)}\n"
        "target_parent_a_fraction\t0.333333\n"
        "target_parent_a_copy_number\t1\n"
        "dstat_parent_a_copy_number\t0.61\n"
        "target_cox1_lineage\tparent_b\n"
        "primary_contig\tpri_chr1\n"
        "reciprocal_contig\trec_chr1\n"
    )


if __name__=="__main__":
    main()
