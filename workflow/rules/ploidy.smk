rule nquire_create:
    input:
        bam=lambda w: f"{RESULTS}/mapping/{PLOIDY_PANEL}/{w.sample}.dedup.bam"
    output:
        bin=f"{RESULTS}/ploidy/nquire/{{sample}}.bin"
    log:
        f"{RESULTS}/logs/nquire_create.{{sample}}.log"
    conda:
        "../../envs/ploidy.yaml"
    params:
        prefix=lambda w: f"{RESULTS}/ploidy/nquire/{w.sample}"
    shell:
        "nQuire create -b {input.bam:q} -o {params.prefix:q} -f 0.15 -c 6 > {log:q} 2>&1"

rule nquire_denoise:
    input:
        bin=rules.nquire_create.output.bin
    output:
        den=f"{RESULTS}/ploidy/nquire/{{sample}}_denoised.bin",
        hist=f"{RESULTS}/ploidy/nquire/{{sample}}.hist.tsv",
        raw_model=f"{RESULTS}/ploidy/nquire/{{sample}}.lrd.raw.tsv",
        den_model=f"{RESULTS}/ploidy/nquire/{{sample}}.lrd.denoised.tsv"
    log:
        f"{RESULTS}/logs/nquire_models.{{sample}}.log"
    conda:
        "../../envs/ploidy.yaml"
    params:
        prefix=lambda w: f"{RESULTS}/ploidy/nquire/{w.sample}"
    shell:
        "(nQuire denoise {input.bin:q} -o {params.prefix}_denoised; "
        "nQuire histo {output.den:q} > {output.hist:q}; "
        "nQuire lrdmodel {input.bin:q} > {output.raw_model:q}; "
        "nQuire lrdmodel {output.den:q} > {output.den_model:q}) 2> {log:q}"

rule ploidy_summary:
    input:
        hists=expand(f"{RESULTS}/ploidy/nquire/{{sample}}.hist.tsv", sample=SAMPLES),
        models=expand(f"{RESULTS}/ploidy/nquire/{{sample}}.lrd.denoised.tsv", sample=SAMPLES),
        samples=rules.validate_inputs.output.normalized
    output:
        summary=f"{RESULTS}/ploidy/ploidy_summary.tsv"
    log:
        f"{RESULTS}/logs/ploidy_summary.log"
    conda:
        "../../envs/core.yaml"
    params:
        histdir=lambda wildcards, input: str(Path(input.hists[0]).parent),
        min_density=config["analysis"].get("ploidy",{}).get("min_spectrum_density",1000),
        quantile=config["analysis"].get("ploidy",{}).get("diploid_ratio_quantile",.95),
        mapping_panel=PLOIDY_PANEL
    shell:
        "python workflow/scripts/summarize_ploidy.py --samples {input.samples:q} "
        "--hist-dir {params.histdir:q} --min-spectrum-density {params.min_density} "
        "--diploid-ratio-quantile {params.quantile} --mapping-panel {params.mapping_panel:q} "
        "--output {output.summary:q} > {log:q} 2>&1"
