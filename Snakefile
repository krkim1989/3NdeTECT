configfile: "config/config.yaml"

include: "workflow/rules/common.smk"
include: "workflow/rules/mapping.smk"
include: "workflow/rules/ploidy.smk"
include: "workflow/rules/mitochondrial.smk"
include: "workflow/rules/ancestry.smk"
include: "workflow/rules/report.smk"

rule all:
    input:
        rules.report.output.html,
        rules.ploidy_summary.output.summary,
        expand(f"{RESULTS}/ancestry/{{panel}}/dosage_summary.tsv", panel=PANELS),
        (rules.reciprocal_concordance.output.tsv if RECIP_REF_KEY else []),
        rules.local_ancestry.output.summary,
        rules.window_ancestry.output.summary,
        rules.dstat.output.tsv,
        rules.dstat.output.calibration,
        rules.dosage_states.output.tsv,
        rules.hybrid_validation.output.tsv,
        (rules.cox1_summary.output if COX1_ENABLED else []),
        (f"{RESULTS}/kmer/kmer_ancestry.tsv" if config["analysis"]["kmer"].get("enabled", False) else []),
