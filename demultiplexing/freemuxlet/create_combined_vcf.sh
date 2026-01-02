
#!/bin/bash

# Combine the VCF files directly
bcftools concat -Oz -o /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/combined.common_rare_phased.snps.vcf.gz /share/hennlab/data/genomes/CAAPA_freeze2_PHASED_common_rare/snps/*common_rare_phased.snps.vcf.gz

bcftools concat -Oz -f /share/hennlab/projects/NCR_scRNAseq/data/common_vcfs/file_list.txt -o /share/hennlab/projects/NCR_scRNAseq/data/common_vcfs/caapa2.vcf.gz

# Index the combined VCF file
bcftools index combined.common_rare_phased.snps.vcf.gz

echo "Combining complete. Output file: combined.common_rare_phased.snps.vcf.gz"


#bcftools view -S /share/hennlab/projects/CAAPA2_functional_annotation/allele_frequencies_and_scores/data/pop_filtered/pop_list/Europ/Nama.csv /share/hennlab/projects/NCR_scRNAseq/data/common_vcfs/caapa2.vcf.gz -Oz -o /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/Nama_only.vcf.gz

# bcftools view -i 'INFO/AF > 0.05' /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/Nama_only.vcf.gz -o /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/Nama_only_5pct.vcf
#thurs -run start from below

#now sort the  file bgzip Nama_only_5pct.vcf 


#downloaded bed from https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2362827405_26b2L1zAes4FzSNRme9OohUWeOHs&clade=mammal&org=&db=hg38&hgta_group=genes&hgta_track=knownGene&hgta_table=knownGene&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=
#bcftools view -R /share/hennlab/projects/NCR_scRNAseq/data/common_vcfs/UCSC_genes_only.bed /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/Nama_only_5pct.vcf -o /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/Nama_5pct_snps_in_genes.vcf
#vcftools --vcf /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/Nama_only_5pct.vcf --bed /share/hennlab/projects/NCR_scRNAseq/data/common_vcfs/UCSC_genes_only.bed --recode --out /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/Nama_5pct_snps_in_genes.vcf
  

  /share/hennlab/projects/NCR_scRNAseq/results/TEStmulti_analysis_GEX2_CSF/outs/per_sample_outs/TEStmulti_analysis_GEX2_CSF/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
  /share/hennlab/projects/NCR_scRNAseq/results/TEStmulti_analysis_GEX_CSF/outs/per_sample_outs/TEStmulti_analysis_GEX_CSF/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz