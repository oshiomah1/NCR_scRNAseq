#!/bin/bash

#SBATCH --partition=production
#SBATCH --job-name=cellbender
#SBATCH --nodes=1 #THIS SHOULD ALWAYS BE One
#SBATCH --ntasks=20 # equivalent to cpus, stick to around 20 max on gc64, or gc128
#SBATCH --mem=250G
#SBATCH --time=4-10:30:00
#SBATCH --output=cellbender%j.out
#SBATCH --error=cellbender%j.err
#SBATCH --mail-user=oyageshio@ucdavis.edu
#SBATCH --mail-type=ALL

# Exit immediately if a command exits with a non-zero status
#set -e

# Load conda environment
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate cellbender

# Define input parameters as variables for reusability
INPUT_DIR="/share/hennlab/projects/NCR_scRNAseq/results/cellranger/batch1/TBSS01/per_sample_outs/TBSS01/count"
INPUT_FILE="${INPUT_DIR}/sample_filtered_feature_bc_matrix.h5"
OUTPUT_DIR="${INPUT_DIR}/cellbender"
OUTPUT_FILE="${OUTPUT_DIR}/TBSS01_filtered.h5"

#EXPECTED_CELLS=500
#TOTAL_DROPLETS=2000

# Ensure the output directory exists
mkdir -p "${OUTPUT_DIR}"

# Run CellBender
cellbender remove-background \
    --input "${INPUT_FILE}" \
    --output "${OUTPUT_FILE}" \
    --fpr  0.1

echo "CellBender processing completed successfully."
