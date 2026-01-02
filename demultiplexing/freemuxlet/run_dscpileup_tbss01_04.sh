#!/bin/bash
#SBATCH --partition=production
#SBATCH --job-name=popscle
#SBATCH --nodes=1 #THIS SHOULD ALWAYS BE One
#SBATCH --ntasks=20 # equivalent to cpus, stick to around 20 max on gc64, or gc128
#SBATCH --mem=250G
#SBATCH --time=4-10:30:00
#SBATCH --output=popscler%j.out
#SBATCH --error=popscle%j.err
#SBATCH --mail-user=oyageshio@ucdavis.edu
#SBATCH --mail-type=ALL

#THIS FOR FOR TBSS01-O4 CSF GEX 

source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate doubletfinder_seurat_jupyter

# Define file paths as variables
#SAM_FILE="/share/hennlab/projects/NCR_scRNAseq/results/TEStmulti_analysis_GEX2_CSF/outs/per_sample_outs/TEStmulti_analysis_GEX2_CSF/count/sample_alignments.bam"
SAM_FILE="/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/filtered_bam/tbssbatch1_filtered.bam"
VCF_FILE="/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence_AF.vcf"
OUTPUT_DIR="/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/dscpileup/tbss01_04_GEX_CSF"
GROUP_LIST="/share/hennlab/projects/NCR_scRNAseq/results/TEStmulti_analysis_GEX2_CSF/outs/per_sample_outs/TEStmulti_analysis_GEX2_CSF/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"

# Run popscle dsc-pileup with variables
popscle dsc-pileup \
    --sam "$SAM_FILE" \
    --tag-group CB \
    --tag-UMI UB \
    --vcf "$VCF_FILE" \
    --out "$OUTPUT_DIR" \
    --group-list "$GROUP_LIST"

#popscle dsc-pileup --sam /share/hennlab/projects/NCR_scRNAseq/results/TEStmulti_analysis_GEX2_CSF/outs/per_sample_outs/TEStmulti_analysis_GEX2_CSF/count/sample_alignments.bam --vcf /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence_AF.vcf --out /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/dscpileup/tbss01_04

#popscle dsc-pileup --sam /share/hennlab/projects/NCR_scRNAseq/results/TEStmulti_analysis_GEX2_CSF/outs/per_sample_outs/TEStmulti_analysis_GEX2_CSF/count/sample_alignments.bam --tag-group CB --tag-UMI UB --vcf /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence_AF.vcf --out /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/dscpileup/tbss01_04_GEX_CSF --group-list /share/hennlab/projects/NCR_scRNAseq/results/TEStmulti_analysis_GEX2_CSF/outs/per_sample_outs/TEStmulti_analysis_GEX2_CSF/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz

 