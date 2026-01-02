#!/bin/bash
#SBATCH --job-name=vireob4wgs #here
#SBATCH --account=genome-center-grp
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH -p high
#SBATCH --mail-user=oyageshio@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/demultiplexing/vireo/script_logs/%x_%A_%a.out
#SBATCH --error=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/demultiplexing/vireo/script_logs/%x_%A_%a.err

set -euo pipefail
module purge
module load conda
eval "$($(command -v mamba) shell hook --shell bash)"
mamba activate /quobyte/bmhenngrp/conda_envs/vireo2

# Variables
BASE_CELL_DIR="/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/demultiplex/cellsnp/batch4" #here
BASE_OUT_DIR="/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/demultiplex/vireob4_with_WGS" #here
WGS_VCF="/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/demultiplex/batch4_donors.sub.vcf.gz" #here
N_EXPECTED_DONORS=26 #here

mkdir -p "$BASE_OUT_DIR"

# Sample IDs
SAMPLES=("TBSS13" "TBSS14" "TBSS15" "TBSS16") #here

for SAMPLE in "${SAMPLES[@]}"; do
    (
        CELL_DIR="${BASE_CELL_DIR}/${SAMPLE}"
        OUT_DIR="${BASE_OUT_DIR}/${SAMPLE}"

        echo "Processing sample ${SAMPLE}"
        echo "Input: ${CELL_DIR}"
        echo "Output: ${OUT_DIR}"

        # Run vireo
        vireo \
          -c "$CELL_DIR" \
          -d "$WGS_VCF" \
          -N "$N_EXPECTED_DONORS" \
          -o "$OUT_DIR"

        echo "✅ Finished ${SAMPLE}"
    ) &
done

wait
echo "All samples processed!"
