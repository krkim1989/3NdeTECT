# Changelog

## 0.2.0 - 2026-08-18

- Validate the effective Snakemake configuration and all reference assets.
- Shard joint genotyping by contig through GenomicsDBImport.
- Calibrate dosage, windows, HMM emissions and k-mer estimates with parental controls.
- Require independent near-fixation in both parental panels for diagnostic SNPs.
- Add beta-binomial local ancestry, constrained/unconstrained model comparison and segment filters.
- Correct D/fd allele orientation, add multiple outgroups, per-individual tests and FDR.
- Add explicit empty-result failures, atomic outputs and software/config provenance.
- Rename the workflow to 3NdeTECT.
- Support different chromosome names for primary and reciprocal assemblies and
  restrict HaplotypeCaller to configured validation intervals.
- Reuse a selected primary/reciprocal BAM for nQuire instead of mapping every
  sample a second time.
- Add a deterministic, lightweight end-to-end dataset with known triploid
  ancestry truth for GitHub smoke testing.
- Fix pandas-version-dependent Boolean filtering in the local-ancestry HMM
  summary.
- Restore the manuscript D-statistic workflow with synthetic copy calibration,
  per-outgroup cross-validation, discrete dosage states, maternal-lineage
  consistency checks and an integrated triploid-hybrid call.
- Add optional FASTQ-to-COX1 extraction with per-lineage independent mapping,
  consensus FASTA, coverage/identity assignment and maternal-lineage integration.
- Extend D-statistic calibration to 0–3 copies, reject non-monotonic/out-of-range
  calibration and propagate its status to integrated output.
- Flag pedigree/COX1 disagreement as `pedigree_cox1_discordant` without treating
  reference-guided COX1 assignment as proof that pedigree metadata is incorrect;
  also validate sample IDs, role/group assignments, selected COX1 samples and
  comparable COX1 reference lengths.
- Quote workflow paths consistently so spaces and shell metacharacters in paths
  cannot break commands.
- Make the distributed configuration and sample-sheet examples taxon-agnostic;
  oyster names remain only as study provenance and optional manuscript guidance.
- Require triploid targets, validate feasible ancestry states and active-method
  thresholds, and make maternal retention checks symmetric for parent A and B.
- Treat the sample-sheet `reference` column as optional provenance and validate
  COX1 maternal keys plus legacy/multiple-outgroup consistency at startup.
- Split nQuire into an isolated environment to avoid its HTSlib conflict, add a
  clean-environment GitHub Actions smoke job, and include Conda >=24.7.1 for
  Snakemake environment inspection.

## 0.1.0 - 2026-08-18

- Initial triploidy and hybrid-ancestry workflow.
