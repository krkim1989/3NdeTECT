# Figures & publication tables.
#
# Every rule here consumes tables the analysis already produced and renders a
# raster figure (PNG for the HTML report + LZW TIFF for publication) plus a
# machine-readable source-data TSV. All labels/ploidy/chromosomes are taken from
# the data or the config, so the same rules serve any cross, not a fixed study.

rule figure_chromosome_painting:
    input:
        segments=rules.local_ancestry.output.segments,
        fai=PRIMARY_FAI
    output:
        figures=fig_outputs("chromosome_painting"),
        source=f"{FIGDIR}/source_data/chromosome_painting.tsv"
    log:
        f"{RESULTS}/logs/figure_chromosome_painting.log"
    conda:
        "../../envs/core.yaml"
    params:
        base=fig_base,
        fmts=fig_formats_arg(),
        label=ANCESTRY_LABEL
    shell:
        "python workflow/scripts/plot_chromosome_painting.py --segments {input.segments:q} "
        "--fai {input.fai:q} --ancestry-label {params.label:q} --focal-copy parent_a "
        "--output-base {params.base:q} --source-data {output.source:q} --formats {params.fmts:q} "
        f"--dpi-png {FIG_DPI_PNG} --dpi-tiff {FIG_DPI_TIFF} > {{log:q}} 2>&1"


rule figure_chromosome_painting_parent_b:
    input:
        segments=rules.local_ancestry.output.segments,
        fai=PRIMARY_FAI
    output:
        figures=fig_outputs("chromosome_painting_parentB"),
        source=f"{FIGDIR}/source_data/chromosome_painting_parentB.tsv"
    log:
        f"{RESULTS}/logs/figure_chromosome_painting_parentB.log"
    conda:
        "../../envs/core.yaml"
    params:
        base=fig_base,
        fmts=fig_formats_arg(),
        label=ANCESTRY_LABEL_B
    shell:
        "python workflow/scripts/plot_chromosome_painting.py --segments {input.segments:q} "
        "--fai {input.fai:q} --ancestry-label {params.label:q} --focal-copy parent_b "
        "--output-base {params.base:q} --source-data {output.source:q} --formats {params.fmts:q} "
        f"--dpi-png {FIG_DPI_PNG} --dpi-tiff {FIG_DPI_TIFF} > {{log:q}} 2>&1"


rule figure_window_ancestry:
    input:
        windows=rules.window_ancestry.output.windows,
        fai=PRIMARY_FAI
    output:
        figures=fig_outputs("window_ancestry"),
        source=f"{FIGDIR}/source_data/window_ancestry.tsv"
    log:
        f"{RESULTS}/logs/figure_window_ancestry.log"
    conda:
        "../../envs/core.yaml"
    params:
        base=fig_base,
        fmts=fig_formats_arg()
    shell:
        "python workflow/scripts/plot_window_ancestry.py --windows {input.windows:q} "
        "--fai {input.fai:q} --output-base {params.base:q} --source-data {output.source:q} "
        f"--formats {{params.fmts:q}} --dpi-png {FIG_DPI_PNG} --dpi-tiff {FIG_DPI_TIFF} > {{log:q}} 2>&1"


rule figure_dstat_calibration:
    input:
        calibration=rules.dstat.output.calibration,
        dstat=rules.dstat.output.tsv
    output:
        figures=fig_outputs("dstat_calibration"),
        source=f"{FIGDIR}/source_data/dstat_calibration.tsv"
    log:
        f"{RESULTS}/logs/figure_dstat_calibration.log"
    conda:
        "../../envs/core.yaml"
    params:
        base=fig_base,
        fmts=fig_formats_arg()
    shell:
        "python workflow/scripts/plot_dstat_calibration.py --calibration {input.calibration:q} "
        "--dstat {input.dstat:q} --output-base {params.base:q} --source-data {output.source:q} "
        f"--formats {{params.fmts:q}} --dpi-png {FIG_DPI_PNG} --dpi-tiff {FIG_DPI_TIFF} > {{log:q}} 2>&1"


rule figure_integer_dosage:
    input:
        states=rules.dosage_states.output.tsv
    output:
        figures=fig_outputs("integer_dosage"),
        source=f"{FIGDIR}/source_data/integer_dosage.tsv"
    log:
        f"{RESULTS}/logs/figure_integer_dosage.log"
    conda:
        "../../envs/core.yaml"
    params:
        base=fig_base,
        fmts=fig_formats_arg(),
        label=ANCESTRY_LABEL
    shell:
        "python workflow/scripts/plot_integer_dosage.py --states {input.states:q} "
        "--ancestry-label {params.label:q} --output-base {params.base:q} --source-data {output.source:q} "
        f"--formats {{params.fmts:q}} --dpi-png {FIG_DPI_PNG} --dpi-tiff {FIG_DPI_TIFF} > {{log:q}} 2>&1"


rule figure_nquire_ploidy:
    input:
        ploidy=rules.ploidy_summary.output.summary
    output:
        figures=fig_outputs("nquire_ploidy"),
        source=f"{FIGDIR}/source_data/nquire_ploidy.tsv"
    log:
        f"{RESULTS}/logs/figure_nquire_ploidy.log"
    conda:
        "../../envs/core.yaml"
    params:
        base=fig_base,
        fmts=fig_formats_arg()
    shell:
        "python workflow/scripts/plot_nquire_ploidy.py --ploidy {input.ploidy:q} "
        "--output-base {params.base:q} --source-data {output.source:q} "
        f"--formats {{params.fmts:q}} --dpi-png {FIG_DPI_PNG} --dpi-tiff {FIG_DPI_TIFF} > {{log:q}} 2>&1"


rule figure_kmer_ancestry:
    input:
        kmer=f"{RESULTS}/kmer/kmer_ancestry.tsv"
    output:
        figures=fig_outputs("kmer_ancestry"),
        source=f"{FIGDIR}/source_data/kmer_ancestry.tsv"
    log:
        f"{RESULTS}/logs/figure_kmer_ancestry.log"
    conda:
        "../../envs/core.yaml"
    params:
        base=fig_base,
        fmts=fig_formats_arg()
    shell:
        "python workflow/scripts/plot_kmer_ancestry.py --kmer {input.kmer:q} "
        "--output-base {params.base:q} --source-data {output.source:q} "
        f"--formats {{params.fmts:q}} --dpi-png {FIG_DPI_PNG} --dpi-tiff {FIG_DPI_TIFF} > {{log:q}} 2>&1"


rule result_tables:
    input:
        dosage=f"{RESULTS}/ancestry/primary/dosage_summary.tsv",
        window=rules.window_ancestry.output.summary,
        dstat=rules.dstat.output.tsv,
        states=rules.dosage_states.output.tsv,
        local=rules.local_ancestry.output.summary,
        ploidy=rules.ploidy_summary.output.summary,
        segments=rules.local_ancestry.output.segments,
        kmer=(f"{RESULTS}/kmer/kmer_ancestry.tsv" if KMER_ENABLED else []),
        reciprocal=(rules.reciprocal_concordance.output.tsv if RECIP_REF_KEY else [])
    output:
        manifest=f"{FIGDIR}/tables/manifest.tsv"
    log:
        f"{RESULTS}/logs/result_tables.log"
    conda:
        "../../envs/core.yaml"
    params:
        outdir=lambda wildcards, output: str(Path(output.manifest).parent),
        kmer=lambda wildcards, input: f"--kmer {shlex.quote(str(input.kmer))}" if input.kmer else "",
        reciprocal=lambda wildcards, input: f"--reciprocal {shlex.quote(str(input.reciprocal))}" if input.reciprocal else ""
    shell:
        "python workflow/scripts/make_tables.py --dosage {input.dosage:q} "
        "--window-summary {input.window:q} --dstat {input.dstat:q} --states {input.states:q} "
        "--local {input.local:q} --ploidy {input.ploidy:q} --segments {input.segments:q} "
        "{params.kmer} {params.reciprocal} --outdir {params.outdir:q} --manifest {output.manifest:q} "
        "> {log:q} 2>&1"


rule figures:
    input:
        figure_targets
    output:
        done=touch(f"{FIGDIR}/figures.done")
