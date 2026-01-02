#!/bin/bash

# Load conda environment
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate cell_snp_lite

# Define constants
VCF_FILE="/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence_AF.vcf"
THREADS=10
MIN_MAF=0.1
MIN_COUNT=20
#make one for batch

# List of sample IDs
#SAMPLES=("TBSS13" "TBSS14" "TBSS15" "TBSS16")
SAMPLES=("TBSS13" "TBSS15" "TBSS16")

# Parallel processing
for SAMPLE in "${SAMPLES[@]}"; do
    (
        # Define paths
        BAM_FILE="/share/hennlab/projects/NCR_scRNAseq/results/cellranger/batch4/$SAMPLE/per_sample_outs/$SAMPLE/count/sample_alignments.bam"
        BARCODES_FILE="/share/hennlab/projects/NCR_scRNAseq/results/cellranger/batch4/$SAMPLE/per_sample_outs/$SAMPLE/count/barcodes.tsv.gz"
        OUTPUT_DIR="/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/cellsnp/batch4/$SAMPLE"
        # Create output directory if it doesn't exist
        mkdir -p $OUTPUT_DIR

        echo "Sample: $SAMPLE"
        echo "BAM_FILE: $BAM_FILE"
        echo "BARCODES_FILE: $BARCODES_FILE"
        echo "OUTPUT_DIR: $OUTPUT_DIR"

        # Run cellsnp-lite
        cellsnp-lite \
            -s "$BAM_FILE" \
            -b "$BARCODES_FILE" \
            -O "$OUTPUT_DIR" \
            -R "$VCF_FILE" \
            -p "$THREADS" \
            --minMAF "$MIN_MAF" \
            --minCOUNT "$MIN_COUNT" \
            --gzip

        echo "$SAMPLE cellsnp complete :)"
    ) &
done

# Wait for all background processes to finish
wait
echo "All samples processed!"
