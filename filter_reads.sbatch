#!/bin/bash
#SBATCH --output=extract_homo_%j.log
#SBATCH --error=extract_homo_%j.err

# Paths
kraken_dir=<OUTPUT_DIR_OF_SBATCH.KRAKEN>
samplesheet=<SAMPLESHEET.CSV>
output_dir=<OUTPUT_DIR_OF_FILTERED_READS>

mkdir -p "$output_dir"

# Skip header line of CSV
tail -n +2 "$samplesheet" | while IFS=',' read -r sample fastq1 fastq2 strandedness; do
    base_kraken="${kraken_dir}/${sample}.kraken"
    if [[ -f "$base_kraken" ]]; then
        echo "Processing $sample..."
        <full_path_to_extract_kraken_reads.py> \
            --report "${kraken_dir}/${sample}.report" \
            -k "${kraken_dir}/${sample}.kraken" \
            -s "$fastq1" \
            -s2 "$fastq2" \
            -o "${output_dir}/${sample}_clean_R1.fq" \
            -o2 "${output_dir}/${sample}_clean_R2.fq" \
            --include-parents --taxid 9606
    else
        echo "⚠️Kraken file not found for sample $sample at $base_kraken"
    fi
done
