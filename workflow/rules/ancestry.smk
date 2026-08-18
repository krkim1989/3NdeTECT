rule diagnostic_sites:
    input:
        vcf=rules.joint_genotype.output.vcf,
        samples=rules.validate_inputs.output.normalized
    output:
        sites=f"{RESULTS}/ancestry/{{panel}}/diagnostic_sites.tsv.gz"
    log:
        f"{RESULTS}/logs/diagnostic_sites.{{panel}}.log"
    conda:
        "../../envs/core.yaml"
    params:
        a=config["analysis"]["parent_a_group"],
        b=config["analysis"]["parent_b_group"],
        strict=config["analysis"]["diagnostic"]["strict_max_minor_af"],
        moderate=config["analysis"]["diagnostic"]["moderate_max_minor_af"],
        min_called=config["analysis"]["diagnostic"]["min_called_alleles"],
        max_missing=config["analysis"]["diagnostic"]["max_missing_fraction"],
        min_dp=config["analysis"]["diagnostic"].get("min_parent_dp", 6),
        min_gq=config["analysis"]["diagnostic"].get("min_parent_gq", 20),
        min_qual=config["analysis"]["diagnostic"].get("min_variant_qual", 30)
    shell:
        "python workflow/scripts/diagnostic_sites.py --vcf {input.vcf:q} --samples {input.samples:q} "
        "--group-a {params.a} --group-b {params.b} --strict {params.strict} "
        "--moderate {params.moderate} --min-called-alleles {params.min_called} "
        "--max-missing-fraction {params.max_missing} --min-parent-dp {params.min_dp} "
        "--min-parent-gq {params.min_gq} --min-qual {params.min_qual} "
        "--output {output.sites:q} > {log:q} 2>&1"

rule ancestry_calls:
    input:
        vcf=rules.joint_genotype.output.vcf,
        sites=rules.diagnostic_sites.output.sites,
        samples=rules.validate_inputs.output.normalized
    output:
        calls=f"{RESULTS}/ancestry/{{panel}}/target_calls.tsv.gz"
    log:
        f"{RESULTS}/logs/ancestry_calls.{{panel}}.log"
    conda:
        "../../envs/core.yaml"
    params:
        target=config["analysis"]["target_group"],
        a=config["analysis"]["parent_a_group"],
        b=config["analysis"]["parent_b_group"],
        min_dp=config["analysis"]["dosage"]["min_dp"],
        min_gq=config["analysis"]["dosage"]["min_gq"]
    shell:
        "python workflow/scripts/extract_ancestry_calls.py --vcf {input.vcf:q} --sites {input.sites:q} "
        "--samples {input.samples:q} --target-group {params.target:q} --group-a {params.a:q} "
        "--group-b {params.b:q} --min-dp {params.min_dp} "
        "--min-gq {params.min_gq} --output {output.calls:q} > {log:q} 2>&1"

rule ancestry_dosage:
    input:
        calls=rules.ancestry_calls.output.calls
    output:
        summary=f"{RESULTS}/ancestry/{{panel}}/dosage_summary.tsv"
    log:
        f"{RESULTS}/logs/ancestry_dosage.{{panel}}.log"
    conda:
        "../../envs/core.yaml"
    shell:
        "python workflow/scripts/summarize_ancestry.py --calls {input.calls:q} --output {output.summary:q} > {log:q} 2>&1"

rule dosage_states:
    input:
        calls=f"{RESULTS}/ancestry/primary/target_calls.tsv.gz",
        cox1=(rules.cox1_summary.output.summary if COX1_ENABLED else [])
    output:
        tsv=f"{RESULTS}/hybrid_validation/dosage_states.tsv"
    log:
        f"{RESULTS}/logs/dosage_states.log"
    conda:
        "../../envs/core.yaml"
    params:
        min_dp=config["analysis"].get("hybrid_validation", {}).get("min_state_dp", 20),
        tolerance=config["analysis"].get("hybrid_validation", {}).get("dosage_state_tolerance", .10),
        max_forbidden=config["analysis"].get("hybrid_validation", {}).get("max_forbidden_fraction", .01),
        maternal=(f"--maternal-parent-group {shlex.quote(config['analysis'].get('hybrid_validation', {}).get('maternal_parent_group'))}" if config["analysis"].get("hybrid_validation", {}).get("maternal_parent_group") else ""),
        parents=f"--parent-a-group {shlex.quote(config['analysis']['parent_a_group'])} --parent-b-group {shlex.quote(config['analysis']['parent_b_group'])}",
        cox1=lambda wildcards, input: f"--cox1-summary {shlex.quote(str(input.cox1))}" if input.cox1 else ""
    shell:
        "python workflow/scripts/dosage_states.py --calls {input.calls:q} --min-dp {params.min_dp} "
        "--tolerance {params.tolerance} --max-forbidden-fraction {params.max_forbidden} {params.parents} {params.maternal} {params.cox1} "
        "--output {output.tsv:q} > {log:q} 2>&1"

rule local_ancestry:
    input:
        calls=f"{RESULTS}/ancestry/primary/target_calls.tsv.gz"
    output:
        segments=f"{RESULTS}/ancestry/local_ancestry_segments.tsv",
        summary=f"{RESULTS}/ancestry/local_ancestry_summary.tsv"
    log:
        f"{RESULTS}/logs/local_ancestry.log"
    conda:
        "../../envs/core.yaml"
    params:
        trans=config["analysis"]["local_ancestry"]["transition_per_mb"],
        minsites=config["analysis"]["local_ancestry"]["min_sites_per_segment"],
        minbp=config["analysis"]["local_ancestry"].get("min_segment_bp", 1000000),
        rho=config["analysis"]["local_ancestry"].get("overdispersion_rho", .02),
        states=",".join(str(x) for x in config["analysis"]["local_ancestry"].get("states", [0,1,2,3])),
        expected=",".join(str(x) for x in config["analysis"]["local_ancestry"].get("expected_fractions", [0,1/3,2/3,1])),
        full="--compare-unconstrained" if config["analysis"]["local_ancestry"].get("compare_unconstrained", True) else ""
    shell:
        "python workflow/scripts/local_ancestry_hmm.py --calls {input.calls:q} "
        "--transition-per-mb {params.trans} --states {params.states:q} --expected-fractions {params.expected:q} "
        "--rho {params.rho} --min-sites {params.minsites} --min-bp {params.minbp} {params.full} "
        "--segments {output.segments:q} --summary {output.summary:q} > {log:q} 2>&1"

rule window_ancestry:
    input:
        calls=f"{RESULTS}/ancestry/primary/target_calls.tsv.gz"
    output:
        windows=f"{RESULTS}/ancestry/window_ancestry.tsv",
        summary=f"{RESULTS}/ancestry/window_ancestry_summary.tsv"
    log:
        f"{RESULTS}/logs/window_ancestry.log"
    conda:
        "../../envs/core.yaml"
    params:
        size=config["analysis"]["windows"]["size"],
        minsites=config["analysis"]["windows"]["min_sites"]
    shell:
        "python workflow/scripts/window_ancestry.py --calls {input.calls:q} --window-size {params.size} "
        "--min-sites {params.minsites} --windows {output.windows:q} --summary {output.summary:q} > {log:q} 2>&1"

rule reciprocal_concordance:
    input:
        primary=f"{RESULTS}/ancestry/primary/dosage_summary.tsv",
        reciprocal=f"{RESULTS}/ancestry/reciprocal/dosage_summary.tsv"
    output:
        tsv=f"{RESULTS}/ancestry/reciprocal_concordance.tsv"
    log:
        f"{RESULTS}/logs/reciprocal_concordance.log"
    conda:
        "../../envs/core.yaml"
    shell:
        "python workflow/scripts/reciprocal_concordance.py --primary {input.primary:q} "
        "--reciprocal {input.reciprocal:q} --output {output.tsv:q} > {log:q} 2>&1"

rule dstat:
    input:
        vcf=f"{RESULTS}/variants/primary/cohort.joint.vcf.gz",
        samples=rules.validate_inputs.output.normalized
    output:
        tsv=f"{RESULTS}/dstat/dstat.tsv",
        calibration=f"{RESULTS}/dstat/synthetic_calibration.tsv"
    log:
        f"{RESULTS}/logs/dstat.log"
    conda:
        "../../envs/core.yaml"
    params:
        a=config["analysis"]["parent_a_group"],
        b=config["analysis"]["parent_b_group"],
        target=config["analysis"]["target_group"],
        out=" ".join(f"--outgroup {shlex.quote(x)}" for x in config["analysis"].get("outgroup_groups", [config["analysis"]["outgroup_group"]])),
        block=config["analysis"]["dstat"]["block_size"],
        min_called=config["analysis"]["dstat"].get("min_called_alleles",4),
        min_dp=config["analysis"]["dstat"].get("min_sample_dp",6),
        min_gq=config["analysis"]["dstat"].get("min_sample_gq",20),
        min_qual=config["analysis"]["dstat"].get("min_variant_qual",30),
        individual="--per-individual" if config["analysis"]["dstat"].get("per_individual",True) else "",
        estimator=config["analysis"]["dstat"].get("estimator", "manuscript_single_orientation"),
        cal_max=config["analysis"]["dstat"].get("calibration_max_copy", 3.0),
        cal_step=config["analysis"]["dstat"].get("calibration_step", .1)
    shell:
        "python workflow/scripts/dstat.py --vcf {input.vcf:q} --samples {input.samples:q} "
        "--p1 {params.b:q} --p2 {params.target:q} --p3 {params.a:q} {params.out} "
        "--block-size {params.block} --min-called-alleles {params.min_called} "
        "--min-sample-dp {params.min_dp} --min-sample-gq {params.min_gq} --min-qual {params.min_qual} "
        "--estimator {params.estimator:q} --calibration-max-copy {params.cal_max} "
        "--calibration-step {params.cal_step} --calibration-output {output.calibration:q} "
        "{params.individual} --output {output.tsv:q} > {log:q} 2>&1"

rule kmer_ancestry:
    input:
        samples=rules.validate_inputs.output.normalized,
        reads1=[row["read1"] for row in SAMPLE_ROWS],
        reads2=[row["read2"] for row in SAMPLE_ROWS]
    output:
        tsv=f"{RESULTS}/kmer/kmer_ancestry.tsv"
    log:
        f"{RESULTS}/logs/kmer_ancestry.log"
    threads: config.get("resources", {}).get("bwa_threads", 8)
    resources:
        mem_mb=config.get("resources", {}).get("mem_mb", 32000)
    conda:
        "../../envs/core.yaml"
    params:
        work=f"{RESULTS}/kmer/work",
        k=config["analysis"]["kmer"].get("k", 21),
        min_count=config["analysis"]["kmer"].get("min_count", 2),
        min_diag=config["analysis"]["kmer"].get("min_diagnostic_kmers",1000),
        memory=config["analysis"]["kmer"].get("memory_gb",2),
        a=config["analysis"]["parent_a_group"],
        b=config["analysis"]["parent_b_group"]
    shell:
        "python workflow/scripts/kmer_ancestry.py --samples {input.samples:q} --work {params.work:q} "
        "--group-a {params.a:q} --group-b {params.b:q} --k {params.k} "
        "--min-count {params.min_count} --min-diagnostic-kmers {params.min_diag} "
        "--threads {threads} --memory {params.memory} --output {output.tsv:q} > {log:q} 2>&1"

rule hybrid_validation:
    input:
        ploidy=rules.ploidy_summary.output.summary,
        dosage=f"{RESULTS}/ancestry/primary/dosage_summary.tsv",
        states=rules.dosage_states.output.tsv,
        local=rules.local_ancestry.output.summary,
        dstat=rules.dstat.output.tsv,
        reciprocal=(rules.reciprocal_concordance.output.tsv if RECIP_REF_KEY else []),
        kmer=(f"{RESULTS}/kmer/kmer_ancestry.tsv" if config["analysis"]["kmer"].get("enabled", False) else []),
        cox1=(rules.cox1_summary.output.summary if COX1_ENABLED else [])
    output:
        tsv=f"{RESULTS}/hybrid_validation/hybrid_validation.tsv"
    log:
        f"{RESULTS}/logs/hybrid_validation.log"
    conda:
        "../../envs/core.yaml"
    params:
        reciprocal=lambda wildcards, input: f"--reciprocal {shlex.quote(str(input.reciprocal))}" if input.reciprocal else "",
        kmer=lambda wildcards, input: f"--kmer {shlex.quote(str(input.kmer))}" if input.kmer else "",
        cox1=lambda wildcards, input: f"--cox1 {shlex.quote(str(input.cox1))}" if input.cox1 else "",
        tier=config["analysis"].get("hybrid_validation", {}).get("tier", "strict"),
        min_copy=config["analysis"].get("hybrid_validation", {}).get("min_parent_a_copy", .10),
        min_methods=config["analysis"].get("hybrid_validation", {}).get("min_methods", 3),
        min_outgroups=config["analysis"].get("hybrid_validation", {}).get("min_significant_outgroups", 1),
        alpha=config["analysis"].get("hybrid_validation", {}).get("alpha", .05),
        expected=config["analysis"].get("hybrid_validation", {}).get("expected_parent_a_copies", 2.0),
        tolerance=config["analysis"].get("hybrid_validation", {}).get("simple_cross_tolerance", .35)
    shell:
        "python workflow/scripts/hybrid_validation.py --ploidy {input.ploidy:q} --dosage {input.dosage:q} "
        "--states {input.states:q} --local {input.local:q} --dstat {input.dstat:q} "
        "{params.reciprocal} {params.kmer} {params.cox1} --tier {params.tier:q} --min-parent-a-copy {params.min_copy} "
        "--min-methods {params.min_methods} --min-significant-outgroups {params.min_outgroups} "
        "--alpha {params.alpha} --expected-parent-a-copies {params.expected} "
        "--simple-cross-tolerance {params.tolerance} --output {output.tsv:q} > {log:q} 2>&1"
