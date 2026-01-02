#!/bin/bash
#SBATCH --job-name=thinvcf5
#SBATCH --account=genome-center-grp
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=12:00:00
#SBATCH -p high
#SBATCH --mail-user=oyageshio@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/demultiplexing/vireo/script_logs/%x_%A_%a.out
#SBATCH --error=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/demultiplexing/vireo/script_logs/%x_%A_%a.err
#SBATCH --mail-user=oyageshio@ucdavis.edu
#SBATCH --mail-type=ALL

#for CELLSNP_VCF since each batch shares a vcf just pick one pwer batch to do it
# paths/quobyte/bmhenngrp/from-lssc0/data/genomes/NCR_Phase_1/final_merged/NCR_Phase1_allchr.filtered.vcf.gz
DONOR_VCF=/quobyte/bmhenngrp/from-lssc0/data/genomes/NCR_Phase_1/final_merged/NCR_Phase1_allchr.filtered.vcf.gz
CELLSNP_VCF=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/demultiplex/cellsnp/batch5/TBSS17/cellSNP.base.vcf.gz #change batch as needed
DONOR_LIST=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/batch5clean.txt #change from b1-b5
OUT_DIR=/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/demultiplex/
mkdir -p $OUT_DIR

module load bcftools

# subset to SNPs seen in scRNAseq (TBSS0x) and only batchx donors
bcftools view $DONOR_VCF \
  -R $CELLSNP_VCF \
  -S $DONOR_LIST \
  -Oz -o $OUT_DIR/batch5_donors.sub.vcf.gz  #change here

module load tabix
# index for vireo
tabix -p vcf $OUT_DIR/batch5_donors.sub.vcf.gz # change here
