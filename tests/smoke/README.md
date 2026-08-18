# End-to-end smoke dataset

This deterministic dataset represents two diploid parent panels, one triploid
target with one parent-A copy and two parent-B copies, and one diploid outgroup.
The two 200-kb references contain the same sequence under different contig names,
which exercises reciprocal-reference interval handling. Expected values are in
`data/truth.tsv`.

The compressed FASTQs and FASTAs are small enough to keep in Git. BWA, FASTA and
GATK indexes are generated locally because they are derived files:

```bash
make smoke-data
make smoke
```

`make smoke` runs the complete workflow and then checks ploidy, dosage,
reciprocal-reference concordance, k-mer ancestry, local ancestry, D-statistic,
FASTQ-derived COX1 consensus/lineage, and report creation against
`data/truth.tsv`.

This dataset tests software integration and estimator direction, not biological
performance on a complex or repetitive genome.
