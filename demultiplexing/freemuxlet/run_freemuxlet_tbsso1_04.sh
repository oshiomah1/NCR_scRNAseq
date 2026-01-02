#!/bin/bash
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate doubletfinder_seurat_jupyter
popscle freemuxlet --plp /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/dscpileup/tbss01_04_GEX_CSF --nsample 29 --out /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/freemuxlet/tbss01_04_GEX_CSF
 