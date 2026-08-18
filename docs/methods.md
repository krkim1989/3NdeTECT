# Methods implemented by 3NdeTECT

## Triploidy validation

Reads are aligned to the reference selected for each biological group. nQuire
creates and denoises the allele-frequency spectrum. The workflow reports density
near 1/3 and 2/3 relative to density near 1/2. Interpretation is made against a
species-matched diploid panel because depth, divergence and heterozygosity affect
the magnitude of the ratio.

## Diagnostic ancestry dosage

Biallelic SNPs are selected from independent diploid `parent_a` and `parent_b`
panels. Target triploids never participate in marker ascertainment. Both panels
must independently pass the near-fixation threshold. Allelic depths are oriented
to parental origin and pooled by sample. Raw ancestry fractions are calibrated
using the observed pure-parent panel endpoints before conversion to copy number.
The reported moderate tier includes the strict markers; the strict tier is its
high-divergence subset.

The complete marker-selection and dosage procedure is repeated on an optional
reciprocal reference. The concordance table reports the per-sample absolute copy
difference rather than silently averaging reference-dependent estimates.

Non-overlapping windows summarize genomic heterogeneity after the same parental
control calibration. The reported variation is descriptive and is not treated
as an independent pedigree estimator.

## Local ancestry

A distance-aware hidden Markov model uses beta-binomial emissions at the
configured parental-A copy fractions. The default constrained states are 0, 1
and 2 copies; an optional unconstrained 0-through-ploidy model is fitted as a
sensitivity analysis. State paths are collapsed into chromosome segments, and
segments below the configured marker-count or base-pair threshold are flagged.
The model describes ancestry tracts and does not by itself identify a unique
backcross generation.

## D statistic

For topology `D(P1, P2, P3, O)`, allele-frequency ABBA and BABA contributions are
summed in physical blocks after per-call depth and genotype-quality filtering.
The manuscript-compatible estimator uses the single outgroup-reference
orientation described in Methods 2.5.6; the optional `frequency_symmetric`
estimator includes the reciprocal outgroup-alternative term. The workflow reports
D and fd with delete-one-block jackknife standard errors, Z scores and
FDR-adjusted P values for pooled and optional per-individual tests.

For every empirical test and outgroup, the same callable sites are reused to
construct synthetic P2 mixtures from P1 and P3 allele frequencies. The configured
P3-copy grid spans the full 0–3 triploid range and is converted to D with the
identical estimator. The observed D is interpolated only when the calibration
curve is monotonic and the observation lies within its range; otherwise the copy
estimate is reported as missing with an explicit calibration status. The value is therefore
dataset/outgroup-specific; an uncalibrated D value is never treated directly as
copy number. Supplying two outgroup groups produces separate calibration curves
and an explicit cross-outgroup validation.

## Integrated triploid-hybrid validation

High-depth diagnostic loci are assigned to the nearest 0, 1, 2 or 3 parental-A
copy state. Clustering near these discrete thirds supports germline ancestry over
diffuse low-level contamination. Maternal retention is evaluated symmetrically:
for a parent-A maternal lineage, zero-copy parent-A tracts are forbidden; for a
parent-B maternal lineage, three-copy parent-A tracts are forbidden. A disagreement
between a configured maternal lineage and an assigned COX1 lineage is reported explicitly
as `pedigree_cox1_discordant` rather than being given a consistency call. This is
a QC flag for disagreement between pedigree metadata and reference-guided COX1
assignment, not proof that the recorded pedigree is incorrect.

The integrated table combines triploidy, diagnostic-SNP dosage, calibrated
D-statistics, local ancestry and optional mapping-independent k-mers. Reciprocal
mapping is retained as a bias sensitivity check. It reports whether hybrid
ancestry is supported by the configured number of independent methods and
whether the observed parental-A copy dosage is compatible with the simple
breeding expectation. It deliberately does not assign a unique backcross
generation without parental broodstock genotypes.

## FASTQ-derived COX1 consensus

When enabled, reads from the configured sample roles are aligned independently
to each supplied COX1 lineage reference. Independent mapping avoids ordering
bias from competing, closely related mitochondrial references. At positions
passing minimum depth and base-quality filters, a majority base is called;
otherwise the consensus contains `N`. Candidate references are ranked by
`breadth × identity`, and the best lineage is reported only when coverage,
identity and score-margin thresholds pass. The selected consensus sequences are
written to FASTA and coverage/identity metrics are retained for inspection.
References must represent comparable COX1 intervals; by default the longest may
not exceed 1.10 times the shortest reference length.

An assigned parent-B COX1 lineage can activate the manuscript maternal-lineage
consistency check, and a parent-A assignment activates the reciprocal check.
COX1 is treated as a
maternal-lineage marker, not as evidence that the entire nuclear genome belongs
to that species. Because this is reference-guided consensus reconstruction rather
than reference-free mitochondrial assembly, discordance should be reviewed for
coverage, reference bias, insufficient lineage divergence and possible NUMTs.

## Mapping-independent k-mers

Meryl builds per-sample canonical spectra after removing singleton errors.
Parent-panel core sets are intersected and subtracted against the opposite
parental pool to obtain conservative diagnostic k-mers. Partial databases are
guarded by completion markers and undersized diagnostic sets fail explicitly.
The optional stage reports raw signals and parent-control-calibrated ancestry.
