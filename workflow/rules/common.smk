import shlex
import sys
from pathlib import Path
from snakemake.utils import validate

sys.path.insert(0, str(Path("workflow").resolve()))
from lib import (
    build_figure_targets,
    discover_contigs,
    figure_outputs,
    load_sample_sheet,
    reference_assets,
    requested_contigs,
    validate_analysis_config,
    validate_contigs,
    validate_cox1_references,
    validate_sample_semantics,
)

validate(config, "../../resources/schema/config.schema.yaml")
validate_analysis_config(config)

SAMPLE_SHEET = config["samples"]
RESULTS = config.get("results", "results")
PRIMARY_REF_KEY = config["analysis"]["primary_reference"]
PRIMARY_REF = config["references"][PRIMARY_REF_KEY]
RECIP_REF_KEY = config["analysis"].get("reciprocal_reference")
PANELS = ["primary"] + (["reciprocal"] if RECIP_REF_KEY else [])
PLOIDY_PANEL = config["analysis"].get("ploidy", {}).get("reference_panel", "primary")
if PLOIDY_PANEL not in PANELS:
    raise ValueError(f"ploidy.reference_panel must be one of {PANELS}: {PLOIDY_PANEL}")
REF_BY_PANEL = {"primary": PRIMARY_REF}
if RECIP_REF_KEY:
    REF_BY_PANEL["reciprocal"] = config["references"][RECIP_REF_KEY]

SAMPLE_ROWS = load_sample_sheet(SAMPLE_SHEET)
validate_sample_semantics(SAMPLE_ROWS, config)

SAMPLES = [r["sample"] for r in SAMPLE_ROWS]
ROW = {r["sample"]: r for r in SAMPLE_ROWS}
TARGETS = [r["sample"] for r in SAMPLE_ROWS if r["role"] == "target"]
PARENT_A = [r["sample"] for r in SAMPLE_ROWS if r["role"] == "parent_a"]
PARENT_B = [r["sample"] for r in SAMPLE_ROWS if r["role"] == "parent_b"]
OUTGROUP = [r["sample"] for r in SAMPLE_ROWS if r["role"] == "outgroup"]

COX1_CONFIG = config["analysis"].get("mitochondrial", {}).get("cox1", {})
COX1_ENABLED = COX1_CONFIG.get("enabled", False)
COX1_REFERENCES = COX1_CONFIG.get("references", {})
COX1_LINEAGES = list(COX1_REFERENCES)
COX1_ROLES = set(COX1_CONFIG.get("sample_roles", ["target"]))
COX1_SAMPLES = [r["sample"] for r in SAMPLE_ROWS if r["role"] in COX1_ROLES]
if COX1_ENABLED:
    if not COX1_SAMPLES:
        raise ValueError(f"COX1 is enabled but sample_roles select no samples: {sorted(COX1_ROLES)}")
    validate_cox1_references(COX1_CONFIG)

REFERENCE_KEY_BY_PANEL = {"primary": PRIMARY_REF_KEY}
if RECIP_REF_KEY:
    REFERENCE_KEY_BY_PANEL["reciprocal"] = RECIP_REF_KEY
REQUESTED_CONTIGS_BY_PANEL = {
    panel: requested_contigs(config["analysis"].get("chromosomes", []), panel, REFERENCE_KEY_BY_PANEL[panel])
    for panel in PANELS
}
CONTIGS_BY_PANEL = {
    panel: (
        validate_contigs(REF_BY_PANEL[panel], REQUESTED_CONTIGS_BY_PANEL[panel])
        if REQUESTED_CONTIGS_BY_PANEL[panel]
        else discover_contigs(REF_BY_PANEL[panel], [])
    )
    for panel in PANELS
}
ALL_REFERENCE_INPUTS = []
for ref in config["references"].values():
    ALL_REFERENCE_INPUTS.extend(reference_assets(ref))

# Shared figure helpers live in the common module so Snakemake's lint rule does
# not see a mixture of helper functions and rules in figures.smk.
FIG_CFG = config.get("report", {}).get("figures", {})
FIGURES_ENABLED = FIG_CFG.get("enabled", True)
FIG_FORMATS = FIG_CFG.get("formats", ["png", "tiff"])
FIG_DPI = FIG_CFG.get("dpi", {}) or {}
FIG_DPI_PNG = FIG_DPI.get("png", 150)
FIG_DPI_TIFF = FIG_DPI.get("tiff", 300)
ANCESTRY_LABEL = FIG_CFG.get("ancestry_label", config["analysis"]["parent_a_group"])
ANCESTRY_LABEL_B = FIG_CFG.get("reciprocal_ancestry_label", config["analysis"]["parent_b_group"])
PAINT_ANCESTRIES = FIG_CFG.get("painting_ancestries", ["parent_a"])
TABLES_ENABLED = FIG_CFG.get("tables", True)
KMER_ENABLED = config["analysis"]["kmer"].get("enabled", False)
FIGDIR = f"{RESULTS}/figures"
PRIMARY_FAI = f"{PRIMARY_REF}.fai"


# Keep Snakefile helpers as lambdas; implementation lives in workflow/lib.py.
fig_outputs = lambda stem: figure_outputs(FIGDIR, FIG_FORMATS, stem)
fig_formats_arg = lambda: ",".join(FIG_FORMATS)
fig_base = lambda wildcards, output: str(Path(output.figures[0]).with_suffix(""))
figure_targets = lambda _: build_figure_targets(
    FIGDIR, FIG_FORMATS, FIGURES_ENABLED, PAINT_ANCESTRIES, KMER_ENABLED, TABLES_ENABLED
)

rule validate_inputs:
    input:
        samples=SAMPLE_SHEET,
        references=ALL_REFERENCE_INPUTS
    output:
        ok=f"{RESULTS}/validation/inputs.ok",
        normalized=f"{RESULTS}/validation/samples.normalized.tsv"
    log:
        f"{RESULTS}/logs/validate_inputs.log"
    params:
        refs=" ".join(
            f"--reference {shlex.quote(key + '=' + value)}"
            for key, value in config["references"].items()
        )
    conda:
        "../../envs/core.yaml"
    shell:
        "python workflow/scripts/validate_inputs.py {params.refs} "
        "--samples {input.samples:q} --output {output.normalized:q} --ok {output.ok:q} > {log:q} 2>&1"
