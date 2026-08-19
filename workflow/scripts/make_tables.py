#!/usr/bin/env python3
"""Curate the analysis outputs into publication-ready result tables (TSV).

Every input is one of the pipeline's existing per-topic TSVs; this step only
selects, renames, and (where useful) aggregates columns into clean tables and
writes them to a single output directory plus a manifest. Optional inputs
(reciprocal, k-mer) are simply skipped when absent, so the same rule serves runs
with or without those analyses. Nothing here is organism-specific.
"""
import argparse
from pathlib import Path

import pandas as pd

from plotting import write_source_data


def read(path):
    return pd.read_csv(path, sep="\t") if path and Path(path).is_file() else None


def keep(df, cols):
    return df[[c for c in cols if c in df.columns]].copy()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dosage", required=True, help="ancestry/primary/dosage_summary.tsv")
    p.add_argument("--window-summary", required=True, help="ancestry/window_ancestry_summary.tsv")
    p.add_argument("--dstat", required=True, help="dstat/dstat.tsv")
    p.add_argument("--states", required=True, help="hybrid_validation/dosage_states.tsv")
    p.add_argument("--local", required=True, help="ancestry/local_ancestry_summary.tsv")
    p.add_argument("--ploidy", required=True, help="ploidy/ploidy_summary.tsv")
    p.add_argument("--segments", help="ancestry/local_ancestry_segments.tsv (per-tract table)")
    p.add_argument("--kmer", help="kmer/kmer_ancestry.tsv (optional)")
    p.add_argument("--reciprocal", help="ancestry/reciprocal_concordance.tsv (optional)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--manifest", required=True)
    a = p.parse_args()

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    produced = []

    def emit(name, df, title):
        if df is None or df.empty:
            return
        path = outdir / name
        write_source_data(df, path)
        produced.append({"table": name, "title": title, "n_rows": len(df)})

    # Table: ploidy evidence
    emit("Table_ploidy_evidence.tsv",
         keep(read(a.ploidy), ["sample", "group", "role", "mapping_panel", "expected_ploidy",
                               "best_model", "thirds_to_half_ratio", "diploid_ratio_threshold",
                               "ploidy_call"]),
         "nQuire ploidy evidence per sample")

    # Table: global ancestry dosage per sample + group aggregate
    dosage = read(a.dosage)
    emit("Table_ancestry_dosage.tsv", dosage, "Genome-wide ancestry dosage per target")
    if dosage is not None and "group" in dosage.columns:
        num = dosage.select_dtypes("number").columns
        agg = dosage.groupby("group")[list(num)].median().reset_index()
        emit("Table_group_ancestry_dosage.tsv", agg, "Median ancestry dosage per group")

    # Table: windowed structure
    emit("Table_window_structure.tsv",
         keep(read(a.window_summary), ["sample", "n_windows", "calibrated_parent_a_fraction",
                                       "sd_across_windows", "expected_binomial_sd",
                                       "overdispersion_ratio", "window_min", "window_max"]),
         "Windowed ancestry heterogeneity")

    # Table: integer dosage states
    states = read(a.states)
    if states is not None:
        frac_cols = [c for c in states.columns if c.startswith("state_") and c.endswith("_fraction")]
        emit("Table_integer_dosage.tsv",
             keep(states, ["sample", "ploidy", "n_sites", "pct_near_integer_dosage",
                           "maternal_lineage_consistency"] + frac_cols),
             "Integer dosage-state composition")

    # Table: D-statistics cross-validation across outgroups
    emit("Table_dstat.tsv",
         keep(read(a.dstat), ["test", "P1", "P2", "P3", "O", "n_sites", "n_blocks",
                              "D", "Z_D", "P_D", "Q_D_FDR", "calibrated_p3_copy_number",
                              "calibration_status"]),
         "D-statistics and calibrated P3 copy number")

    # Table: local-ancestry summary
    emit("Table_local_ancestry.tsv",
         keep(read(a.local), ["sample", "ploidy", "n_segments", "n_pass_segments",
                              "marker_span_bp", "bpweighted_parent_a_ancestry",
                              "filtered_bpweighted_parent_a_ancestry",
                              "unconstrained_3copy_bp_fraction"]),
         "Local-ancestry segment summary")

    # Table: per-individual tracts (constrained segments)
    seg = read(a.segments) if a.segments else None
    if seg is not None and "model" in seg.columns:
        seg = seg[seg.model == "constrained"]
    emit("Table_per_individual_tracts.tsv",
         keep(seg, ["sample", "CHROM", "start", "end", "bp_length", "A_COPY",
                    "n_sites", "mean_A_read_fraction", "pass_filter"]) if seg is not None else None,
         "Per-individual ancestry tracts")

    # Optional tables
    emit("Table_kmer_ancestry.tsv",
         keep(read(a.kmer), ["sample", "group", "role", "ploidy", "raw_parent_a_fraction",
                             "calibrated_parent_a_fraction", "calibrated_parent_a_copy_number"])
         if a.kmer else None,
         "Mapping-independent k-mer ancestry")
    emit("Table_reciprocal_concordance.tsv", read(a.reciprocal) if a.reciprocal else None,
         "Reciprocal-reference concordance")

    write_source_data(pd.DataFrame(produced), a.manifest)
    print(f"wrote {len(produced)} tables to {outdir}")


if __name__ == "__main__":
    main()
