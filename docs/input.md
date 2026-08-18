# Input preparation

3NdeTECT uses file-based configuration. Put reference and result paths plus
analysis options in `config/config.yaml`, and put one sample per row with FASTQ
paths in `config/samples.tsv`. The command line only selects the YAML and the
available cores; FASTQ paths do not need to be repeated on the command line.
All taxon names are user-defined. `parent_a`, `parent_b`, `target` and `outgroup`
are analytical roles rather than oyster-specific biological labels.
The default triploid local-ancestry state space is 0, 1, 2 and 3 parent-A
copies. Restrict it only when an independently justified biological model
excludes a state.

- Use paired, gzip-compressed FASTQ files. Sample IDs may contain only letters,
  numbers, dot, underscore and hyphen.
- Use chromosome-level references when chromosome painting is intended.
- Every reference requires pre-built BWA (`.bwt`), FASTA (`.fai`) and GATK
  sequence-dictionary (`.dict`) indexes. Keeping reference preparation outside
  the workflow avoids concurrent writes to shared reference directories.
- Include several unrelated diploid individuals per parental panel. Ten per
  species were used in the oyster study.
- Every `target` row must have `ploidy: 3`. Non-target controls may use another
  positive ploidy, while both parental diagnostic-marker panels must be diploid.
- The outgroup must be sufficiently diverged to orient allele sharing and should
  be tested for sensitivity when a second outgroup is available.
- To reproduce the oyster manuscript cross-validation, add both `HK` and `CV`
  outgroup rows to the sample sheet and set `outgroup_groups: [HK, CV]`. Set
  `hybrid_validation.min_significant_outgroups: 2` to require agreement.
- Set `hybrid_validation.maternal_parent_group` only when mitochondrial or
  pedigree evidence identifies the maternal lineage. Use `null` when unknown;
  the maternal nuclear-consistency field will then be reported as not evaluated.
  If pedigree metadata and the reference-guided COX1 assignment disagree, the
  workflow reports the QC state `pedigree_cox1_discordant`; it does not declare
  either source incorrect automatically.
- For FASTQ-derived COX1 assignment, provide at least two single-record COX1
  FASTA files under `analysis.mitochondrial.cox1.references`. Reference keys
  should match the corresponding sample-sheet group names when the result will
  be used for maternal parent-A/parent-B consistency checks. All references
  should cover the same homologous interval; `max_reference_length_ratio`
  controls the fail-fast length check (default 1.10).
- Keep `dstat.calibration_max_copy: 3.0` for general triploid validation. Values
  outside a monotonic calibration curve are reported as missing rather than
  clamped to an endpoint.
- The sample-sheet `reference` column is optional provenance only. It does not
  choose a per-sample mapper reference; primary and reciprocal references apply
  to the complete cohort.

When primary and reciprocal assemblies use different contig names, configure
`analysis.chromosomes` as a mapping with separate `primary` and `reciprocal`
lists. A single list remains supported when names are shared.
