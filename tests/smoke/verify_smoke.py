#!/usr/bin/env python3
"""Verify that the synthetic end-to-end run recovers its known truth."""

import argparse
import csv
from pathlib import Path


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def one(path, **wanted):
    matches = [row for row in rows(path) if all(row.get(key) == value for key, value in wanted.items())]
    if len(matches) != 1:
        raise AssertionError(f"expected one matching row in {path}, found {len(matches)}: {wanted}")
    return matches[0]


def near(label, observed, expected, tolerance=0.05):
    if abs(float(observed) - float(expected)) > tolerance:
        raise AssertionError(f"{label}: observed {observed}, expected {expected} +/- {tolerance}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default="tests/smoke/results")
    parser.add_argument("--truth", default="tests/smoke/data/truth.tsv")
    args = parser.parse_args()
    result = Path(args.results)
    truth = {row["quantity"]: row["value"] for row in rows(Path(args.truth))}
    expected_copy = float(truth["target_parent_a_copy_number"])
    expected_fraction = float(truth["target_parent_a_fraction"])

    ploidy = one(result / "ploidy/ploidy_summary.tsv", sample="T3N01")
    if ploidy["ploidy_call"] != "triploid_supported":
        raise AssertionError(f"ploidy call was {ploidy['ploidy_call']!r}")

    primary = one(result / "ancestry/primary/dosage_summary.tsv", sample="T3N01", tier="strict")
    reciprocal = one(result / "ancestry/reciprocal/dosage_summary.tsv", sample="T3N01", tier="strict")
    if int(primary["n_sites"]) < 100:
        raise AssertionError("fewer than 100 diagnostic markers survived")
    near("primary dosage copy", primary["calibrated_parent_a_copy_number"], expected_copy)
    near("reciprocal dosage copy", reciprocal["calibrated_parent_a_copy_number"], expected_copy)

    concordance = one(result / "ancestry/reciprocal_concordance.tsv", sample="T3N01", tier="strict")
    if float(concordance["absolute_copy_difference"]) > 0.02:
        raise AssertionError("primary/reciprocal copy estimates are discordant")

    kmer = one(result / "kmer/kmer_ancestry.tsv", sample="T3N01")
    near("k-mer copy", kmer["calibrated_parent_a_copy_number"], expected_copy)

    local = one(result / "ancestry/local_ancestry_summary.tsv", sample="T3N01")
    near("local ancestry", local["bpweighted_parent_a_ancestry"], expected_fraction)
    windows = one(result / "ancestry/window_ancestry_summary.tsv", sample="T3N01")
    near("window ancestry", windows["calibrated_parent_a_fraction"], expected_fraction)

    dstat = one(result / "dstat/dstat.tsv", test="T3N_pooled")
    if int(dstat["n_blocks"]) < 3 or float(dstat["D"]) <= 0:
        raise AssertionError("D-statistic did not recover positive parent-A sharing")
    near("D-statistic calibrated copy", dstat["calibrated_p3_copy_number"], truth["dstat_parent_a_copy_number"], .10)
    if dstat["calibration_status"] != "within_range":
        raise AssertionError(f"D-statistic calibration status was {dstat['calibration_status']!r}")
    hybrid = one(result / "hybrid_validation/hybrid_validation.tsv", sample="T3N01")
    if hybrid["hybrid_call"] != "triploid_hybrid_supported":
        raise AssertionError(f"integrated hybrid call was {hybrid['hybrid_call']!r}")
    if hybrid["simple_cross_compatible"] != "False":
        raise AssertionError("1A:2B smoke target was incorrectly accepted as the configured 2A:1B cross")
    if hybrid["cox1_matches_configured_maternal_group"] != "True":
        raise AssertionError("COX1 did not match the configured maternal group")
    if hybrid["maternal_group_evaluated"] != "parent_b" or float(hybrid["forbidden_opposite_parent_state"]) != 3:
        raise AssertionError("parent-B maternal retention did not evaluate the three-copy parent-A state")
    cox1 = one(result / "cox1/cox1_summary.tsv", sample="T3N01")
    if cox1["cox1_status"] != "assigned" or cox1["cox1_lineage"] != truth["target_cox1_lineage"]:
        raise AssertionError(f"COX1 lineage was not recovered: {cox1}")
    if not (result / "cox1/cox1_consensus.fa").is_file():
        raise AssertionError("selected COX1 consensus FASTA is missing")
    if not (result / "report/index.html").is_file():
        raise AssertionError("HTML report is missing")

    print(
        "PASS: triploid_supported; "
        f"dosage={float(primary['calibrated_parent_a_copy_number']):.4f} copies; "
        f"kmer={float(kmer['calibrated_parent_a_copy_number']):.4f} copies; "
        f"local={float(local['bpweighted_parent_a_ancestry']):.4f}; "
        f"D={float(dstat['D']):.4f}; "
        f"D-copy={float(dstat['calibrated_p3_copy_number']):.4f}; "
        f"COX1={cox1['cox1_lineage']}"
    )


if __name__ == "__main__":
    main()
