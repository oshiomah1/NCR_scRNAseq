# scode to reorder vcf to bam file

 
#module load picard-tools
 
 
 #start from scartch again and see f this works

 # step 1 - use GATK to UPDATESEQUENCE DICTIONARY

gatk UpdateVCFSequenceDictionary -V /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/Nama_5pct_snps_in_genes._v2vcf.recode.vcf --source-dictionary /share/hennlab/projects/NCR_scRNAseq/data/common_vcfs/TBSS01.bam  --output /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.vcf --replace
#grep -n "chr19" /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.vcf | tail -n 1 | cut -d: -f1

 #sed -n '5078857p' /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.vcf 
# step 2 then USE picadools to sort
picard SortVcf I=/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.vcf O=/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence.vcf 
 #grep -n "chr19" /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence.vcf | tail -n 1 | cut -d: -f1

 #sed -n '3118397p' /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence.vcf  

 #step 3 annotate the vcf with AF field

bcftools +fill-tags /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence.vcf -Ov -o /share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/new/Nama_5pct_snps_in_genes.sorted.header.sorted.sequence_AF.vcf  

 