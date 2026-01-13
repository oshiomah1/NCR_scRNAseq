#!/bin/bash
#SBATCH --job-name=multi_db
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --array=1-4
#SBATCH --output=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/logs/%x_%j.out
#SBATCH --error=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/logs/%x_%j.err
#SBATCH --mail-user=oshiomah1234@gmail.com
#SBATCH --mail-type=ALL

module load conda
eval "$($(command -v mamba) shell hook --shell bash)"
mamba activate /quobyte/bmhenngrp/conda_envs/seuratt

set -euo pipefail

Rscript /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/multi_db_run.R \
  /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/raw_merge_all_batches_harm_annotated_all_res12_pca15_noann.rds \
  /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/ScType_multiDB_out_res12_pca15
