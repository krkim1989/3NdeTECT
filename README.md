# 3NdeTECT

Taxon-agnostic workflow for (1) validating triploidy from whole-genome reads and
(2) estimating interspecific ancestry in triploid organisms. It was developed
from analyses of *Crassostrea gigas*–*C. angulata* triploid oysters, but no
species name is hard-coded: parental, target, outgroup and COX1-lineage labels
are supplied entirely by the user.

## What it does

1. validates the sample sheet and reference files;
2. maps reads to primary and optional reciprocal references and calls mixed-ploidy gVCFs/VCFs;
3. reuses the selected ancestry-mapping BAM to test ploidy with nQuire spectra;
4. builds species-diagnostic SNPs **only from diploid reference panels**;
5. estimates triploid ancestry dosage from allele depths;
6. reconstructs configured parental-copy tracts with a distance-aware beta-binomial HMM;
7. tests excess allele sharing using block-jackknifed D and fd statistics;
8. optionally confirms ancestry with mapping-independent diagnostic 21-mers;
9. optionally extracts a COX1 consensus directly from FASTQ reads and assigns
   the maternal mitochondrial lineage against user-supplied references;
10. writes machine-readable QC tables and an HTML report containing optional
   reciprocal-reference and k-mer checks.

The workflow deliberately reports ancestry and pedigree separately. It can show
that a triploid has mosaic ancestry, but it does not label a specific backcross
generation unless parental broodstock or a separately validated demographic
model supports that conclusion.

## Quick start

Requirements: Linux, Conda/Mamba and Snakemake >= 8.

```bash
git clone https://github.com/USERNAME/3ndetect.git
cd 3ndetect
micromamba create -y -n 3ndetect-core -f envs/core.yaml
micromamba activate 3ndetect-core
cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
# edit both files, then validate before launching expensive jobs
snakemake --configfile config/config.yaml --lint
snakemake --configfile config/config.yaml --use-conda --cores 24 --dry-run
snakemake --configfile config/config.yaml --use-conda --cores 24
```

Developer checks can be run with `make check`. The repository CI executes the
same synthetic ancestry and HMM regression tests on every push and pull request.
Run `make smoke` for a real BWA/GATK/nQuire/Meryl end-to-end execution on the
fully synthetic dataset under `tests/smoke`. nQuire uses the separate `envs/ploidy.yaml`
environment because its older HTSlib dependency conflicts with the current
samtools/bcftools versions in `envs/core.yaml`. Snakemake creates both rule
environments automatically with `--use-conda`.

For a SLURM cluster:

```bash
snakemake --configfile config/config.yaml --use-conda \
  --workflow-profile profiles/slurm
```

## Sample sheet

Required columns:

| column | meaning |
|---|---|
| `sample` | unique ID using only letters, numbers, dot, underscore or hyphen |
| `group` | population/species group |
| `role` | `target`, `parent_a`, `parent_b`, `outgroup`, or `control` |
| `ploidy` | target must be `3`; diploid parental marker panels must be `2` |
| `read1`, `read2` | paired FASTQ paths |

An optional `reference` column may be retained as informational provenance, but
it does not select the mapping reference. Mapping references are selected by
`analysis.primary_reference` and `analysis.reciprocal_reference` for the whole
cohort.

At least two diploid parental panels (`parent_a`, `parent_b`) are required for
diagnostic ancestry markers. Target samples must not be used to ascertain those
markers.

## Main outputs

```text
results/
  validation/                 input and software checks
  mapping/                    BAM and mapping QC
  ploidy/                     nQuire spectra and ploidy summary
  variants/                   per-sample gVCFs and joint VCF
  ancestry/
    primary/diagnostic_sites.tsv.gz
    primary/target_calls.tsv.gz
    primary/dosage_summary.tsv
    reciprocal/dosage_summary.tsv
    reciprocal_concordance.tsv
    window_ancestry_summary.tsv
    local_ancestry_segments.tsv
    local_ancestry_summary.tsv
  dstat/dstat.tsv
  dstat/synthetic_calibration.tsv
  hybrid_validation/dosage_states.tsv
  hybrid_validation/hybrid_validation.tsv
  cox1/cox1_summary.tsv       FASTQ-derived lineage, coverage and identity
  cox1/cox1_consensus.fa      selected per-sample consensus sequences
  kmer/kmer_ancestry.tsv      optional
  report/index.html
```

## Interpretation guardrails

- Triploidy is classified from allele-balance modes near 1/3 and 2/3 and the
  nQuire likelihood fit, calibrated against species/reference-matched diploid
  controls. Low-information samples are reported as `insufficient_data`.
- Diagnostic SNPs must be selected from independent diploid panels.
- D statistics test allele-sharing asymmetry. 3NdeTECT generates a
  dataset/outgroup-specific synthetic mixture curve before expressing D on a
  parental-copy scale.
- Mitochondrial lineage identifies the maternal lineage, not a completely pure
  maternal nuclear genome.
- HMM state transitions are ancestry boundaries, not automatically individual
  meiotic recombination events.

See [docs/methods.md](docs/methods.md) and [docs/input.md](docs/input.md).

## Reproducibility and citation

Conda environments pin the major bioinformatics tools. Each result directory
contains logs, and the final report records the exact effective configuration
and software versions. If this workflow is used in a publication, cite the
individual tools listed in `CITATION.cff` and the associated biological study
when available.

## License

MIT for workflow code. Reference assemblies and sequencing data retain their
original licenses and are not redistributed. All versioned test reads,
genotypes, and references in this repository are synthetic.
