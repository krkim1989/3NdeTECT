# 3NdeTECT

3NdeTECT is a taxon-agnostic Snakemake workflow for validating triploidy from
whole-genome sequencing and estimating interspecific ancestry in triploid
organisms. It was developed from *Crassostrea gigas*–*C. angulata* oyster
analyses, but species names, parental groups, target groups, outgroups, and
mitochondrial lineages are supplied by the user rather than hard-coded.

The workflow answers two separate questions:

1. Is the target genuinely triploid?
2. If so, how much of each parental ancestry is present, and is the pattern
   compatible with a simple cross or with older introgression/backcrossing?

Ancestry evidence does not by itself prove a specific pedigree generation.
3NdeTECT reports the genomic pattern and whether it is compatible with the
configured simple-cross expectation; broodstock genotypes or an independently
validated demographic model are needed to label an exact backcross generation.

## Analysis modules

| module | purpose | principal evidence |
|---|---|---|
| Input validation | fail early on inconsistent metadata or missing assets | sample roles, ploidy, FASTQ and index checks |
| Mixed-ploidy calling | call diploid panels and triploid targets together | ploidy-aware GATK gVCF genotypes |
| Triploidy | distinguish diploid-like and triploid-like allele balance | nQuire spectra and diploid-control calibration |
| Diagnostic dosage | estimate parental copies in each target | allele depth at independently ascertained parental SNPs |
| Local ancestry | find ancestry tracts and mosaic structure | distance-aware beta-binomial HMM |
| D and fd | test excess allele sharing | block jackknife, per-individual tests, multiple outgroups |
| Reciprocal mapping | quantify reference-mapping sensitivity | copy-estimate concordance between assemblies |
| k-mer ancestry | optional mapping-independent validation | parental diagnostic k-mers |
| COX1 | infer the maternal mitochondrial lineage | independent mapping to user-supplied COX1 candidates |
| Integrated call | combine independent evidence without hiding conflicts | machine-readable QC table and HTML report |

Diagnostic markers are selected only from diploid parental panels. Target
hybrids are never used to ascertain their own ancestry markers.

## How input is supplied

3NdeTECT is file-configured:

- `config/config.yaml` contains reference paths, group names, thresholds,
  enabled modules, resources, and the result directory.
- `config/samples.tsv` contains one row per sample and the paired FASTQ paths.
- The command line selects the YAML file and compute resources. FASTQ paths do
  not need to be repeated in every command.

This makes a run reproducible and avoids long commands containing dozens of
sample paths.

## Requirements

- Linux
- Conda or Micromamba
- Snakemake 8.30 (provided by `envs/core.yaml`)
- Conda 24.7.1 or newer for Snakemake's rule-environment handling
- paired-end FASTQ data
- one primary reference assembly and, optionally, a reciprocal reference

The core environment contains BWA, samtools, bcftools, GATK, Meryl, Python, and
Snakemake. nQuire is isolated in `envs/ploidy.yaml` because its older HTSlib
dependency conflicts with the current samtools/bcftools stack.

## Installation

Micromamba is recommended:

```bash
git clone https://github.com/krkim1989/3NdeTECT.git
cd 3NdeTECT
micromamba create -y -n 3ndetect-core -f envs/core.yaml
micromamba activate 3ndetect-core
```

Conda can be used instead:

```bash
conda env create -n 3ndetect-core -f envs/core.yaml
conda activate 3ndetect-core
```

Confirm the installation before configuring real data:

```bash
snakemake --version
conda --version
make check
```

`make check` performs Python compilation, unit tests, Snakemake lint, and a
dry-run. `make smoke` additionally runs the complete 70-job workflow on fully
synthetic reads and known ancestry truth. No real-study reads or genotypes are
distributed in this repository.

## Prepare reference assemblies

Each nuclear reference requires BWA indexes, a FASTA index, and a GATK sequence
dictionary. Build them once before starting the workflow:

```bash
bwa index /data/reference/parent_b.fa
samtools faidx /data/reference/parent_b.fa
gatk CreateSequenceDictionary \
  -R /data/reference/parent_b.fa \
  -O /data/reference/parent_b.dict
```

Repeat for the reciprocal reference if used. The `.dict` basename must match
the FASTA basename. Keep reference indexing outside Snakemake so concurrent
jobs cannot attempt to modify the same shared reference.

For COX1, provide at least two single-record FASTA files covering the same
homologous interval. They may represent any candidate maternal lineages; the
workflow is not restricted to *C. angulata*.

## Configure samples

Create working copies of the examples:

```bash
cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
```

Required sample-sheet columns:

| column | meaning |
|---|---|
| `sample` | unique ID containing letters, numbers, `.`, `_`, or `-` |
| `group` | user-defined population or species label |
| `role` | `target`, `parent_a`, `parent_b`, `outgroup`, or `control` |
| `ploidy` | target must be `3`; diagnostic parental panels must be `2` |
| `read1`, `read2` | paired FASTQ paths |
| `reference` | optional provenance label; it does not select the mapper reference |

Example:

```text
sample  group      role      ploidy  read1                    read2
T3N01   Hybrid3N   target    3       /data/T3N01_R1.fastq.gz  /data/T3N01_R2.fastq.gz
PA01    SpeciesA   parent_a  2       /data/PA01_R1.fastq.gz   /data/PA01_R2.fastq.gz
PA02    SpeciesA   parent_a  2       /data/PA02_R1.fastq.gz   /data/PA02_R2.fastq.gz
PB01    SpeciesB   parent_b  2       /data/PB01_R1.fastq.gz   /data/PB01_R2.fastq.gz
PB02    SpeciesB   parent_b  2       /data/PB02_R1.fastq.gz   /data/PB02_R2.fastq.gz
O01     Outgroup1  outgroup  2       /data/O01_R1.fastq.gz    /data/O01_R2.fastq.gz
```

Use several unrelated diploid individuals per parental panel. Small panels can
run technically but make diagnostic-marker ascertainment and calibration less
robust.

## Configure the analysis

Important `config/config.yaml` fields:

| key | function |
|---|---|
| `samples` | path to the sample sheet |
| `results` | output directory |
| `references` | named nuclear-reference FASTA paths |
| `analysis.primary_reference` | reference key used for the main analysis |
| `analysis.reciprocal_reference` | optional second reference key; use `null` to disable |
| `analysis.parent_a_group`, `analysis.parent_b_group` | sample-sheet groups used as parental panels |
| `analysis.target_group` | triploid group |
| `analysis.outgroup_groups` | one or more outgroup group names |
| `analysis.chromosomes` | contigs to process; empty means all contigs |
| `analysis.hybrid_validation.maternal_parent_group` | pedigree/mtDNA maternal prior, or `null` |
| `analysis.mitochondrial.cox1.references` | lineage label to COX1 FASTA mapping |
| `analysis.kmer.enabled` | enable or disable mapping-independent k-mer validation |
| `resources` | per-rule threads and memory requests |

Group labels must match the sample sheet exactly, including capitalization.
When primary and reciprocal assemblies use different contig names, configure
separate `primary` and `reciprocal` chromosome lists.

## Recommended execution order

Always validate and dry-run before launching expensive jobs:

```bash
# 1. Static workflow checks
snakemake --configfile config/config.yaml --lint

# 2. Resolve the DAG without running jobs
snakemake --configfile config/config.yaml \
  --use-conda --cores 1 --dry-run

# 3. Optionally create all Conda environments first
snakemake --configfile config/config.yaml \
  --use-conda --cores 1 --conda-create-envs-only

# 4. Run
snakemake --configfile config/config.yaml \
  --use-conda --cores 24 --printshellcmds --show-failed-logs
```

Snakemake resumes from completed outputs. Re-running the same command does not
normally repeat successful jobs.

For SLURM:

```bash
snakemake --configfile config/config.yaml --use-conda \
  --workflow-profile profiles/slurm
```

## Main outputs and interpretation

| output | interpretation |
|---|---|
| `validation/inputs.ok` | configuration and input assets passed validation |
| `ploidy/ploidy_summary.tsv` | triploidy support or insufficient information |
| `ancestry/primary/dosage_summary.tsv` | calibrated parent-A copy estimate |
| `ancestry/reciprocal_concordance.tsv` | sensitivity to mapping reference |
| `ancestry/window_ancestry_summary.tsv` | ancestry fraction by genomic window |
| `ancestry/local_ancestry_segments.tsv` | HMM ancestry segments and copy states |
| `dstat/dstat.tsv` | ABBA/BABA, D, fd, jackknife uncertainty, and FDR |
| `dstat/synthetic_calibration.tsv` | dataset-specific D-to-copy calibration curve |
| `hybrid_validation/dosage_states.tsv` | discrete triploid parental-copy states |
| `hybrid_validation/hybrid_validation.tsv` | integrated evidence and conflict flags |
| `cox1/cox1_summary.tsv` | maternal-lineage assignment, breadth, identity, margin |
| `kmer/kmer_ancestry.tsv` | optional mapping-independent copy estimate |
| `report/index.html` | human-readable summary report |

Key interpretation rules:

- `triploid_supported` requires informative allele-balance spectra relative to
  matched diploid controls. `insufficient_data` is not evidence of diploidy.
- A parental-copy estimate is meaningful only when enough independent
  diagnostic SNPs pass QC and reciprocal-reference estimates are concordant.
- D statistics test allele-sharing asymmetry. Copy numbers derived from D use a
  dataset/outgroup-specific synthetic calibration and are not universal units.
- `calibration_status=nonmonotonic`, `below_range`, or `above_range` should not
  be converted into a pedigree statement.
- COX1 identifies a mitochondrial lineage, not a completely pure maternal
  nuclear genome.
- HMM boundaries indicate ancestry-state changes, not automatically individual
  meiotic recombination events or an exact number of generations.
- Agreement among dosage, D, k-mer, reciprocal mapping, and tract structure is
  stronger than any single method.

## Troubleshooting and recovery

Start with the failed rule name and its log:

```bash
snakemake --configfile config/config.yaml \
  --use-conda --cores 1 --dry-run

find results/logs -maxdepth 1 -type f -print
tail -n 100 results/logs/<failed-rule>.log
```

After fixing the cause, resume safely:

```bash
snakemake --configfile config/config.yaml \
  --use-conda --cores 24 --rerun-incomplete \
  --printshellcmds --show-failed-logs
```

Common failures:

| symptom | likely cause | action |
|---|---|---|
| `LockException` | an earlier Snakemake process ended unexpectedly | confirm no Snakemake process is running, then run `snakemake --unlock --configfile config/config.yaml` |
| `MissingInputException` | wrong FASTQ/reference path or missing index | check absolute paths and build `.bwt`, `.fai`, and `.dict` files |
| `Conda must be version 24.7.1 or later` | an old base Conda is being used | activate `3ndetect-core` and confirm `which conda` and `conda --version` |
| Conda solve/download failure | transient network, channel, or cache problem | retry environment creation; confirm access to conda-forge, bioconda, and HCC |
| `nQuire: command not found` | running installed-tools mode without the ploidy environment | use normal `--use-conda`; `smoke-installed` requires the ploidy environment on `PATH` |
| reference dictionary error | FASTA and `.dict` basenames or contigs disagree | rebuild the dictionary from the exact FASTA configured in YAML |
| `contig ... not found` | chromosome names differ between assemblies/config | inspect FASTA `.fai` files and use separate primary/reciprocal contig lists |
| BAM/VCF truncated or incomplete | a job was interrupted or storage filled | check free space and `samtools quickcheck`; resume with `--rerun-incomplete` |
| Java out-of-memory / exit 137 | memory request is too small or host killed the job | increase `resources.mem_mb`, reduce concurrent jobs, and rerun |
| `fewer than three informative blocks` | too few callable D-stat sites or blocks | verify group labels/outgroups, relax only justified QC filters, or analyze more sequence |
| D calibration `nonmonotonic` | the chosen data/outgroup do not give an ordered synthetic curve | do not force a copy estimate; inspect outgroup choice, missingness, and allele orientation |
| COX1 `unassigned` | low breadth/identity for all candidates | check read quality and provide a closer, homologous COX1 reference |
| COX1 `ambiguous` | top candidates have nearly equal scores | add discriminating references or sequence; do not force a maternal call |
| `pedigree_cox1_discordant` | configured maternal prior and COX1 disagree | verify metadata and references; report the conflict rather than overriding either source |
| ploidy `insufficient_data` | too few informative heterozygous sites | check depth/mapping and matched diploid controls; do not interpret as diploid |
| too few diagnostic SNPs | small/divergent parental panels, missing data, or overly strict thresholds | verify parental roles first, then review depth/GQ/missingness thresholds |
| primary/reciprocal discordance | reference bias or nonhomologous contig selection | inspect mapping QC and ensure equivalent genomic scope in both assemblies |

Do not delete the entire result directory as a first response. Snakemake can
usually reuse valid outputs. Remove or regenerate only an output proven to be
corrupt, and keep the rule log for diagnosis.

If a run still fails, collect the following before opening an issue:

```bash
snakemake --version
conda --version
git rev-parse HEAD
snakemake --configfile config/config.yaml --summary
```

Include the failed rule, its log, the relevant YAML section with private paths
redacted, and whether the failure reproduces with `--cores 1`.

## Reproducibility and citation

Each result directory contains logs, the effective configuration snapshot, and
software versions. Keep these with the machine-readable tables when archiving
an analysis. Cite the tools listed in `CITATION.cff` and the associated
biological study when applicable.

See [docs/input.md](docs/input.md) for input-design details and
[docs/methods.md](docs/methods.md) for the statistical methods.

## License

MIT for workflow code. Reference assemblies and sequencing data retain their
original licenses and are not redistributed. All versioned test reads,
genotypes, and references in this repository are synthetic.
