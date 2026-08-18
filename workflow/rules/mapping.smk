rule map_reads:
    input:
        ok=rules.validate_inputs.output.ok,
        ref=lambda w: REF_BY_PANEL[w.panel],
        fai=lambda w: REF_BY_PANEL[w.panel] + ".fai",
        bwt=lambda w: REF_BY_PANEL[w.panel] + ".bwt",
        dct=lambda w: str(Path(REF_BY_PANEL[w.panel]).with_suffix("")) + ".dict",
        r1=lambda w: ROW[w.sample]["read1"],
        r2=lambda w: ROW[w.sample]["read2"]
    output:
        bam=f"{RESULTS}/mapping/{{panel}}/{{sample}}.dedup.bam",
        bai=f"{RESULTS}/mapping/{{panel}}/{{sample}}.dedup.bam.bai",
        metrics=f"{RESULTS}/mapping/{{panel}}/{{sample}}.markdup_metrics.txt"
    log:
        f"{RESULTS}/logs/map.{{panel}}.{{sample}}.log"
    threads: config.get("resources", {}).get("bwa_threads", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", 16000)
    conda:
        "../../envs/core.yaml"
    params:
        sort_threads=config.get("resources", {}).get("sort_threads", 4),
        rg=lambda w: f"@RG\\tID:{w.sample}\\tSM:{w.sample}\\tPL:ILLUMINA"
    shell:
        "(bwa mem -t {threads} -R '{params.rg}' {input.ref:q} {input.r1:q} {input.r2:q} | "
        "samtools sort -@ {params.sort_threads} -o {output.bam:q}.tmp -; "
        "gatk MarkDuplicates -I {output.bam:q}.tmp -O {output.bam:q}.markdup -M {output.metrics:q} "
        "--CREATE_INDEX false; mv {output.bam:q}.markdup {output.bam:q}; "
        "samtools index -@ {params.sort_threads} {output.bam:q}; samtools quickcheck -v {output.bam:q}; "
        "rm -f {output.bam:q}.tmp) > {log:q} 2>&1"

rule haplotypecaller:
    input:
        ref=lambda w: REF_BY_PANEL[w.panel],
        bam=rules.map_reads.output.bam,
        bai=rules.map_reads.output.bai
    output:
        vcf=f"{RESULTS}/variants/{{panel}}/gvcf/{{sample}}.g.vcf.gz",
        tbi=f"{RESULTS}/variants/{{panel}}/gvcf/{{sample}}.g.vcf.gz.tbi"
    log:
        f"{RESULTS}/logs/haplotypecaller.{{panel}}.{{sample}}.log"
    threads: config.get("resources", {}).get("gatk_threads", 4)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", 16000)
    conda:
        "../../envs/core.yaml"
    params:
        ploidy=lambda w: int(ROW[w.sample]["ploidy"]),
        intervals=lambda w: " ".join(f"-L {shlex.quote(contig)}" for contig in REQUESTED_CONTIGS_BY_PANEL[w.panel])
    shell:
        "gatk HaplotypeCaller -R {input.ref:q} -I {input.bam:q} -O {output.vcf:q} "
        "-ERC GVCF --sample-ploidy {params.ploidy} {params.intervals} "
        "--native-pair-hmm-threads {threads} > {log:q} 2>&1"

rule genomicsdb_sample_map:
    input:
        vcfs=lambda w: expand(f"{RESULTS}/variants/{w.panel}/gvcf/{{sample}}.g.vcf.gz", sample=SAMPLES),
        tbis=lambda w: expand(f"{RESULTS}/variants/{w.panel}/gvcf/{{sample}}.g.vcf.gz.tbi", sample=SAMPLES)
    output:
        tsv=f"{RESULTS}/variants/{{panel}}/sample_map.tsv"
    log:
        f"{RESULTS}/logs/sample_map.{{panel}}.log"
    conda:
        "../../envs/core.yaml"
    params:
        samples=",".join(SAMPLES)
    shell:
        "python workflow/scripts/make_sample_map.py --samples {params.samples:q} "
        "--vcfs {input.vcfs:q} --output {output.tsv:q} > {log:q} 2>&1"

rule genomicsdb_import:
    input:
        ref=lambda w: REF_BY_PANEL[w.panel],
        sample_map=rules.genomicsdb_sample_map.output.tsv
    output:
        db=directory(f"{RESULTS}/variants/{{panel}}/genomicsdb/{{contig}}")
    log:
        f"{RESULTS}/logs/genomicsdb.{{panel}}.{{contig}}.log"
    threads: config.get("resources", {}).get("gatk_threads", 4)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", 32000)
    conda:
        "../../envs/core.yaml"
    shell:
        "gatk GenomicsDBImport --sample-name-map {input.sample_map:q} "
        "--genomicsdb-workspace-path {output.db:q} -L {wildcards.contig:q} "
        "--reader-threads {threads} > {log:q} 2>&1"

rule genotype_contig:
    input:
        ref=lambda w: REF_BY_PANEL[w.panel],
        db=rules.genomicsdb_import.output.db
    output:
        vcf=f"{RESULTS}/variants/{{panel}}/joint/{{contig}}.vcf.gz",
        tbi=f"{RESULTS}/variants/{{panel}}/joint/{{contig}}.vcf.gz.tbi"
    log:
        f"{RESULTS}/logs/genotype.{{panel}}.{{contig}}.log"
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", 32000)
    conda:
        "../../envs/core.yaml"
    params:
        db_uri=lambda w, input: f"gendb://{input.db}"
    shell:
        "gatk GenotypeGVCFs -R {input.ref:q} -V {params.db_uri:q} "
        "-L {wildcards.contig:q} -O {output.vcf:q} > {log:q} 2>&1"

rule joint_genotype:
    input:
        vcfs=lambda w: expand(
            f"{RESULTS}/variants/{w.panel}/joint/{{contig}}.vcf.gz",
            contig=CONTIGS_BY_PANEL[w.panel],
        ),
        tbis=lambda w: expand(
            f"{RESULTS}/variants/{w.panel}/joint/{{contig}}.vcf.gz.tbi",
            contig=CONTIGS_BY_PANEL[w.panel],
        )
    output:
        vcf=f"{RESULTS}/variants/{{panel}}/cohort.joint.vcf.gz",
        tbi=f"{RESULTS}/variants/{{panel}}/cohort.joint.vcf.gz.tbi"
    log:
        f"{RESULTS}/logs/joint_concat.{{panel}}.log"
    threads: config.get("resources", {}).get("sort_threads", 4)
    conda:
        "../../envs/core.yaml"
    shell:
        "bcftools concat --threads {threads} -a -Oz -o {output.vcf:q} {input.vcfs:q} "
        "> {log:q} 2>&1; bcftools index --threads {threads} -t {output.vcf:q} >> {log:q} 2>&1"
