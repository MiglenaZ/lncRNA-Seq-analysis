#!/bin/bash
# Convert FASTA-like reads (with '>' headers, no quality lines) into valid FASTQ files (with '@' headers and fake quality lines).

INPUT_DIR=<OUTPUT_DIR_OF_FILTERED_READS>
OUTPUT_DIR="${INPUT_DIR}/fixed_fastq"

mkdir -p "$OUTPUT_DIR"

echo "Converting FASTA-like files in $INPUT_DIR to FASTQ.gz format..."

for f in ${INPUT_DIR}/*_clean_R*.fq.gz; do
    [ -e "$f" ] || continue  # Skip if no match
    base=$(basename "$f" .fq.gz)
    out="${OUTPUT_DIR}/${base}_fixed.fq.gz"

    echo "Processing $base ..."
    zcat "$f" | \
    awk 'BEGIN{OFS="\n"}
         /^>/ {sub(/^>/,"@"); h=$0; getline; seq=$0;
               print h,seq,"+",
               substr("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",1,length(seq))
              }' | gzip > "$out"

    if [[ $? -eq 0 ]]; then
        echo "✅ Converted: $out"
    else
        echo "❌ Error converting: $f"
    fi
done
echo "All done! Converted files are in: $OUTPUT_DIR"
