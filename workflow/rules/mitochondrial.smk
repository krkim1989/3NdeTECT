rule prepare_cox1_reference:
    input:
        source=lambda w: COX1_REFERENCES[w.lineage]
    output:
        fa=f"{RESULTS}/cox1/references/{{lineage}}.fa",
        amb=f"{RESULTS}/cox1/references/{{lineage}}.fa.amb",
        ann=f"{RESULTS}/cox1/references/{{lineage}}.fa.ann",
        bwt=f"{RESULTS}/cox1/references/{{lineage}}.fa.bwt",
        pac=f"{RESULTS}/cox1/references/{{lineage}}.fa.pac",
        sa=f"{RESULTS}/cox1/references/{{lineage}}.fa.sa",
        fai=f"{RESULTS}/cox1/references/{{lineage}}.fa.fai"
    log:
        f"{RESULTS}/logs/cox1_reference.{{lineage}}.log"
    conda:
        "../../envs/core.yaml"
    shell:
        "(python workflow/scripts/normalize_cox1_reference.py --input {input.source:q} "
        "--lineage {wildcards.lineage:q} --output {output.fa:q}; "
        "bwa index {output.fa:q}; samtools faidx {output.fa:q}) > {log:q} 2>&1"

rule map_cox1_reads:
    input:
        ref=rules.prepare_cox1_reference.output.fa,
        bwt=rules.prepare_cox1_reference.output.bwt,
        r1=lambda w: ROW[w.sample]["read1"],
        r2=lambda w: ROW[w.sample]["read2"]
    output:
        bam=f"{RESULTS}/cox1/mapping/{{sample}}.{{lineage}}.bam",
        bai=f"{RESULTS}/cox1/mapping/{{sample}}.{{lineage}}.bam.bai"
    log:
        f"{RESULTS}/logs/cox1_map.{{sample}}.{{lineage}}.log"
    threads: config.get("resources", {}).get("bwa_threads", 8)
    conda:
        "../../envs/core.yaml"
    params:
        sort_threads=config.get("resources", {}).get("sort_threads", 4),
        rg=lambda w: f"@RG\\tID:{w.sample}.cox1.{w.lineage}\\tSM:{w.sample}\\tPL:ILLUMINA"
    shell:
        "(bwa mem -t {threads} -R '{params.rg}' {input.ref:q} {input.r1:q} {input.r2:q} | "
        "samtools sort -@ {params.sort_threads} -o {output.bam:q}; "
        "samtools index -@ {params.sort_threads} {output.bam:q}; samtools quickcheck -v {output.bam:q}) > {log:q} 2>&1"

rule cox1_candidate_consensus:
    input:
        bam=rules.map_cox1_reads.output.bam,
        bai=rules.map_cox1_reads.output.bai,
        ref=rules.prepare_cox1_reference.output.fa,
        fai=rules.prepare_cox1_reference.output.fai
    output:
        metrics=f"{RESULTS}/cox1/candidates/{{sample}}.{{lineage}}.metrics.tsv",
        fasta=f"{RESULTS}/cox1/candidates/{{sample}}.{{lineage}}.consensus.fa"
    log:
        f"{RESULTS}/logs/cox1_consensus.{{sample}}.{{lineage}}.log"
    conda:
        "../../envs/core.yaml"
    params:
        min_depth=COX1_CONFIG.get("min_depth", 3),
        min_bq=COX1_CONFIG.get("min_base_quality", 20),
        min_major=COX1_CONFIG.get("min_major_fraction", .60)
    shell:
        "python workflow/scripts/cox1_consensus.py --bam {input.bam:q} --reference {input.ref:q} "
        "--sample {wildcards.sample:q} --lineage {wildcards.lineage:q} --min-depth {params.min_depth} "
        "--min-base-quality {params.min_bq} --min-major-fraction {params.min_major} "
        "--metrics {output.metrics:q} --fasta {output.fasta:q} > {log:q} 2>&1"

rule cox1_summary:
    input:
        metrics=expand(f"{RESULTS}/cox1/candidates/{{sample}}.{{lineage}}.metrics.tsv",sample=COX1_SAMPLES,lineage=COX1_LINEAGES),
        fasta=expand(f"{RESULTS}/cox1/candidates/{{sample}}.{{lineage}}.consensus.fa",sample=COX1_SAMPLES,lineage=COX1_LINEAGES)
    output:
        summary=f"{RESULTS}/cox1/cox1_summary.tsv",
        consensus=f"{RESULTS}/cox1/cox1_consensus.fa"
    log:
        f"{RESULTS}/logs/cox1_summary.log"
    conda:
        "../../envs/core.yaml"
    params:
        min_breadth=COX1_CONFIG.get("min_breadth", .80),
        min_identity=COX1_CONFIG.get("min_identity", .95),
        min_margin=COX1_CONFIG.get("min_score_margin", .001)
    shell:
        "python workflow/scripts/summarize_cox1.py --metrics {input.metrics:q} "
        "--min-breadth {params.min_breadth} --min-identity {params.min_identity} "
        "--min-margin {params.min_margin} --summary {output.summary:q} "
        "--consensus {output.consensus:q} > {log:q} 2>&1"
