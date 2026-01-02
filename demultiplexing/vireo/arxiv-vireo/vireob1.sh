#!/bin/bash

# Load conda environment
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate scRNAseq
#vireosnp                  0.5.8            pyh7cba7a3_0    bioconda
# Define variables
BASE_CELL_DIR="/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/cellsnp/batch1"
BASE_OUT_DIR="/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/vireo_batch_1"
n_donor=29


# List of sample IDs to process
############ edit all your sample IDs here#######################
SAMPLES=("TBSS01" "TBSS02" "TBSS03" "TBSS04") #  edit all your sample IDs here
################################################################################


# Loop through each sample in parallel
for SAMPLE in "${SAMPLES[@]}"; do
    (
        CELL_DIR="${BASE_CELL_DIR}/${SAMPLE}"
        OUT_DIR="${BASE_OUT_DIR}/${SAMPLE}"

        # Create output directory if it doesn't exist
        #mkdir -p "$OUT_DIR"

        
        echo " $SAMPLE"
        echo " $CELL_DIR"
        echo " $OUT_DIR"
        

        # Run vireo
        vireo -c "$CELL_DIR" -N "$n_donor" -o "$OUT_DIR"
        echo "Processed ${SAMPLE}"
    ) &
done

# Wait for all background jobs to finish
wait
echo "Processed all samples"