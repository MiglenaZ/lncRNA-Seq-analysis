#!/bin/bash
#SBATCH --output=rnaseq_%j.log
#SBATCH --error=rnaseq_%j.err

source ~/.bashrc
conda activate env_nf

# gtf downloaded gencode.v49.annotation.gtf.gz
# fasta downloaded GRCh38.p14.genome.fa.gz

nextflow run nf-core/lncpipe -r dev \
        --input <PATH_TO_SAMPLESHEET_AFTER_FIX_FASTA_TO_FASTQ.SH> \
        --outdir <OUTPUT_DIR> \
        --gtf <FULL_PATH_TO_GTF_FILE> \
        --fasta <FULL_PATH_TO_FASTA_FILE> \
        -profile <APPTAINER/DOCKER> \
        -with-mpi
