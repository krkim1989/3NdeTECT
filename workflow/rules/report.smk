rule config_snapshot:
    output:
        yaml=f"{RESULTS}/provenance/config.used.yaml"
    log:
        f"{RESULTS}/logs/config_snapshot.log"
    conda:
        "../../envs/core.yaml"
    script:
        "../scripts/snapshot_config.py"

rule software_versions:
    output:
        tsv=f"{RESULTS}/provenance/software_versions.tsv"
    log:
        f"{RESULTS}/logs/software_versions.log"
    conda:
        "../../envs/core.yaml"
    shell:
        "python workflow/scripts/collect_versions.py --output {output.tsv:q} > {log:q} 2>&1"

rule report:
    input:
        ploidy=rules.ploidy_summary.output.summary,
        dosage=f"{RESULTS}/ancestry/primary/dosage_summary.tsv",
        local=rules.local_ancestry.output.summary,
        windows=rules.window_ancestry.output.summary,
        dstat=rules.dstat.output.tsv,
        dcal=rules.dstat.output.calibration,
        states=rules.dosage_states.output.tsv,
        hybrid=rules.hybrid_validation.output.tsv,
        reciprocal=(rules.reciprocal_concordance.output.tsv if RECIP_REF_KEY else []),
        kmer=(f"{RESULTS}/kmer/kmer_ancestry.tsv" if config["analysis"]["kmer"].get("enabled", False) else []),
        cox1=(rules.cox1_summary.output.summary if COX1_ENABLED else []),
        config=rules.config_snapshot.output.yaml,
        versions=rules.software_versions.output.tsv,
        figures=(rules.figures.output.done if FIGURES_ENABLED else [])
    output:
        html=f"{RESULTS}/report/index.html"
    log:
        f"{RESULTS}/logs/report.log"
    conda:
        "../../envs/core.yaml"
    params:
        title=config.get("report", {}).get("title", "3NdeTECT report"),
        reciprocal=lambda wildcards, input: f"--reciprocal {shlex.quote(str(input.reciprocal))}" if input.reciprocal else "",
        kmer=lambda wildcards, input: f"--kmer {shlex.quote(str(input.kmer))}" if input.kmer else "",
        cox1=lambda wildcards, input: f"--cox1 {shlex.quote(str(input.cox1))}" if input.cox1 else "",
        figures=(f"--figures-dir {shlex.quote(FIGDIR)}" if FIGURES_ENABLED else "")
    shell:
        "python workflow/scripts/make_report.py --hybrid {input.hybrid:q} --ploidy {input.ploidy:q} "
        "--dosage {input.dosage:q} --states {input.states:q} --local {input.local:q} "
        "--windows {input.windows:q} --dstat {input.dstat:q} --dcal {input.dcal:q} "
        "--config {input.config:q} --versions {input.versions:q} "
        "{params.reciprocal} {params.kmer} {params.cox1} {params.figures} "
        "--title {params.title:q} --output {output.html:q} > {log:q} 2>&1"
