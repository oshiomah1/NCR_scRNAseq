#!/bin/bash
#SBATCH --job-name=res1.0pcs15
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --output=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/logs/%x_%j.out
#SBATCH --error=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/logs/%x_%j.err
#SBATCH --mail-user=oshiomah1234@gmail.com
#SBATCH --mail-type=ALL

# --- Safety & environment ---
set -euo pipefail

# Create log dir if missing
mkdir -p /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/logs


# # Use node-local scratch for temp files (faster; optional)
export TMPDIR=${TMPDIR:-/scratch/${USER}/${SLURM_JOB_ID}}
mkdir -p "$TMPDIR"

#use mamba env that has Seurat installed

module load conda
eval "$($(command -v mamba) shell hook --shell bash)"
mamba activate /quobyte/bmhenngrp/conda_envs/seuratt

# HARD LIMIT THREADING TO SLURM REQUEST
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK


# # --- Paths ---
R_SCRIPT="/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/step_3_HARMONY.R"
#R_SCRIPT="/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/step3a_tune-params.R"
#R_SCRIPT="/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/step_3_manual_sctype.R"


# echo "[$(date)] Starting step_3_HARMONY.R on ${SLURM_NODELIST}"
# echo "JobID: ${SLURM_JOB_ID}  CPUs: ${SLURM_CPUS_PER_TASK}  Mem: ${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-unknown}}"

# --- Run ---
Rscript --vanilla "$R_SCRIPT"

echo "[$(date)] Done."


 
