#!/bin/bash

# prebuilt database (https://benlangmead.github.io/aws-indexes/k2) Standard-16 .tar.gz file in July, 2025 collection (https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20250714.tar.gz) was downloaded and unzipped into <DATABASE_DIRECTORY>, but any prebuilt or manually built Kraken2 databases work.
DB_PATH=<DATABASE_DIRECTORY>
OUTDIR=<OUTPUT_DIRECTORY>
SAMPLESHEET=<SAMPLESHEET>

mkdir -p "$OUTDIR"
sed -i 's/\r$//' "$SAMPLESHEET"

# skip header
tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample fastq1 fastq2 strandedness
do
    echo "🔍 Processing $sample ..."
    echo "FASTQ1: $fastq1"
    echo "FASTQ2: $fastq2"

    kraken2 \
        --db "$DB_PATH" \
        --paired "$fastq1" "$fastq2" \
        --threads 8 \
        --use-names \
        --gzip-compressed \
        --memory-mapping \
        --report "$OUTDIR/${sample}.report" \
        --output "$OUTDIR/${sample}.kraken"

    echo "✅ Finished $sample"
done
