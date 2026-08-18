#!/usr/bin/env bash
set -euo pipefail

python tests/smoke/generate_smoke_data.py --output tests/smoke/data
for reference in tests/smoke/data/primary.fa tests/smoke/data/reciprocal.fa; do
    stem=${reference%.fa}
    rm -f "$reference.amb" "$reference.ann" "$reference.bwt" \
      "$reference.pac" "$reference.sa" "$reference.fai" "$stem.dict"
    bwa index "$reference"
    samtools faidx "$reference"
    gatk CreateSequenceDictionary -R "$reference" -O "$stem.dict"
done
