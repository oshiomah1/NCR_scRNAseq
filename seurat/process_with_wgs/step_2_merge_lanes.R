library(Seurat)

library(ggplot2)
library(tibble)
library(purrr)
library(stringr)

batch_paths_list <- list(
  batch1_paths = c(
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch1/TBSS01/TBSS01_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch1/TBSS02/TBSS02_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch1/TBSS03/TBSS03_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch1/TBSS04/TBSS04_qc_only_vdj.rds"
  ),
  batch2_paths = c(
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch2/TBSS05/TBSS05_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch2/TBSS06/TBSS06_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch2/TBSS07/TBSS07_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch2/TBSS08/TBSS08_qc_only_vdj.rds"
  ),
  batch3_paths = c(
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch3/TBSS09/TBSS09_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch3/TBSS10/TBSS10_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch3/TBSS11/TBSS11_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch3/TBSS12/TBSS12_qc_only_vdj.rds"
  ),
  batch4_paths = c(
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch4/TBSS13/TBSS13_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch4/TBSS14/TBSS14_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch4/TBSS15/TBSS15_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch4/TBSS16/TBSS16_qc_only_vdj.rds"
  ),
  batch5_paths = c(
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch5/TBSS17/TBSS17_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch5/TBSS18/TBSS18_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch5/TBSS19/TBSS19_qc_only_vdj.rds",
    "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch5/TBSS20/TBSS20_qc_only_vdj.rds"
  )
)

# batch_paths_list <- list(
#   batch4_paths = c(
#     "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch4/TBSS13/TBSS13_qc_only_vdj.rds",
#     "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch4/TBSS14/TBSS14_qc_only_vdj.rds",
#     "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch4/TBSS15/TBSS15_qc_only_vdj.rds",
#     "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch4/TBSS16/TBSS16_qc_only_vdj.rds"
#   ),
#   batch5_paths = c(
#     "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch5/TBSS17/TBSS17_qc_only_vdj.rds",
#     "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch5/TBSS18/TBSS18_qc_only_vdj.rds",
#     "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch5/TBSS19/TBSS19_qc_only_vdj.rds",
#     "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_wgs/batch5/TBSS20/TBSS20_qc_only_vdj.rds"
#   )
# )


outdir <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

process_batch <- function(batch_paths, batch_name) {
  log_path <- file.path(outdir, paste0(batch_name, "_output.txt"))
  pdf_path <- file.path(outdir, paste0(batch_name, "_plots.pdf"))
  rds_path <- file.path(outdir, paste0(batch_name, "_combined_wdj.rds"))
  
  sink(log_path, split = TRUE)
  on.exit({ try(sink(), silent = TRUE) }, add = TRUE)
  
  pdf(pdf_path)
  on.exit({ if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE) }, add = TRUE)
  
  message("Reading RDS for ", batch_name, " ...")
  seurat_list <- lapply(batch_paths, readRDS)
  
  message("Merging ...")
  combined <- Reduce(function(x, y) merge(x, y), seurat_list)
  combined$batch <- batch_name
  
  ncells <- ncol(combined)
  print(paste("Number of cells in", batch_name, "is", ncells))
  
  print(VlnPlot(combined, features = "nFeature_RNA", group.by = "orig.ident") + guides(fill = "none"))
  print(VlnPlot(combined, features = "nCount_RNA",   group.by = "orig.ident") + guides(fill = "none"))
  if ("percent.mt"   %in% colnames(combined@meta.data)) print(VlnPlot(combined, features = "percent.mt",   group.by = "orig.ident") + guides(fill = "none"))
  if ("nCount_ADT"   %in% colnames(combined@meta.data)) print(VlnPlot(combined, features = "nCount_ADT",   group.by = "orig.ident") + guides(fill = "none"))
  if ("nFeature_ADT" %in% colnames(combined@meta.data)) print(VlnPlot(combined, features = "nFeature_ADT", group.by = "orig.ident") + guides(fill = "none"))
  if (all(c("nCount_RNA","percent.mt") %in% colnames(combined@meta.data))) print(FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt"))
  if (all(c("nCount_ADT","nFeature_ADT") %in% colnames(combined@meta.data))) print(FeatureScatter(combined, feature1 = "nCount_ADT", feature2 = "nFeature_ADT"))
  if ("percent.ribo" %in% colnames(combined@meta.data)) print(VlnPlot(combined, features = "percent.ribo", group.by = "orig.ident") + guides(fill = "none"))
  
  message("Saving ...")
  saveRDS(combined, rds_path)
  
  rm(combined, seurat_list); gc()
  tibble(batch = batch_name, cells = ncells, rds = rds_path, pdf = pdf_path, log = log_path)
}

# run batches one-by-one
summaries <- imap_dfr(batch_paths_list, process_batch)
summaries
