#!/bin/bash

popscle dsc-pileup --sam /share/hennlab/projects/NCR_scRNAseq/results/TEStmulti_analysis_GEX_CSF/outs/per_sample_outs/TEStmulti_analysis_GEX_CSF/count/sample_alignments.bam --vcf /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence_AF.vcf --out /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/dscpileup/test


# script below is to run freemuxlet after popscle
 
#popscle freemuxlet --plp /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/dscpileup/test --out /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/freemuxlet/test --nsample 29