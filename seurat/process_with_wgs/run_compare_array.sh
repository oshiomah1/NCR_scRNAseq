#!/bin/bash
#SBATCH --job-name=compare_res
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --array=1-4
#SBATCH --output=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/logs/%x_%j.out
#SBATCH --error=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/logs/%x_%j.err
#SBATCH --mail-user=oshiomah1234@gmail.com
#SBATCH --mail-type=ALL



module load conda
eval "$($(command -v mamba) shell hook --shell bash)"
mamba activate /quobyte/bmhenngrp/conda_envs/seuratt

set -euo pipefail
 

RDS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/rds_list.txt)
OUTDIR="/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/compare_resolutions"
mkdir -p "$OUTDIR"

 
Rscript /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/process_one_rds.R "$RDS" "$OUTDIR"
