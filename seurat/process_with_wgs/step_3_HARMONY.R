pkgs <- c("Seurat","harmony","ggplot2","dplyr","scales","readxl","stringr","tibble", "openxlsx")
print(sapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "MISSING")))
invisible(lapply(pkgs, require, character.only=TRUE))
cat("All required libraries loaded.\n")
#library(openxlsx)
library(HGNChelper)
library(future)
library(clustree)

plan("sequential")  # force all future-based code to run in the main R session
options(future.globals.maxSize = Inf)  # optional, but harmless if plan is sequential

log_file   <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/log_file_rawmerge_allres12_pca35.txt"
output_pdf <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/raw_merge_output_file_allres12_pca35.pdf"
save_rds   <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/raw_merge_all_batches_harm_annotated_all_res12_pca35_noann.rds"

sink(log_file, split = TRUE); on.exit({try(sink(), TRUE)}, add=TRUE)
pdf(output_pdf, width=8, height=6); on.exit({if (dev.cur()>1) try(dev.off(), TRUE)}, add=TRUE)

paths <- c(
"/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch1_paths_combined_wdj.rds",
"/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch2_paths_combined_wdj.rds",
"/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch3_paths_combined_wdj.rds",
"/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch4_paths_combined_wdj.rds",
"/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch5_paths_combined_wdj.rds")

# 1) Load & merge
seurat_list <- lapply(paths, readRDS)
merged_obj  <- Reduce(function(x,y) merge(x,y), seurat_list)
rm(seurat_list); gc()
cat("Cells:", ncol(merged_obj), "\n")
cat("Batches:", paste(unique(merged_obj$batch), collapse=", "), "\n")
DefaultAssay(merged_obj) <- "RNA"
 
 

# 2) RNA preprocessing

merged_obj <- NormalizeData(merged_obj, normalization.method="LogNormalize", verbose=FALSE) %>%
  FindVariableFeatures(selection.method="vst", nfeatures=2000, verbose=FALSE) %>%
  ScaleData(verbose=FALSE) %>%
  RunPCA(npcs=50, reduction.name="pca", reduction.key="PCA_", verbose=FALSE)

print(ElbowPlot(merged_obj, ndims = 60))

# Before Harmony
p1 <- DimPlot(merged_obj, reduction = "pca", group.by = "batch")

# 3) Harmony on batch
merged_obj <- RunHarmony(
  object = merged_obj,
  group.by.vars = "batch",
  dims.use = 1:35,
  assay.use = "RNA",
  plot_convergence = TRUE
)


# 4) Graph/UMAP on Harmony
merged_obj <- FindNeighbors(merged_obj, reduction="harmony", dims=1:35, verbose=FALSE) %>%
  FindClusters(resolution=1.2, verbose=FALSE) %>%
  RunUMAP(reduction="harmony", dims=1:35,
          reduction.name="rna.umap", reduction.key="rnaUMAP_", verbose=FALSE)


# After Harmony
p2 <- DimPlot(merged_obj, reduction = "harmony", group.by = "batch")
p1
p2


# 5) (Optional) ADT for marker support (not used for Harmony)
if ("ADT" %in% Assays(merged_obj)) {

  merged_obj <- NormalizeData(merged_obj, assay="ADT",
                              normalization.method="CLR", margin=2)
  merged_obj <- ScaleData(merged_obj, assay="ADT")

  VariableFeatures(merged_obj[["ADT"]]) <- rownames(merged_obj[["ADT"]])

  merged_obj <- RunPCA(merged_obj, assay="ADT",
                       reduction.name="apca", reduction.key="APCA_")

  merged_obj <- RunUMAP(merged_obj, reduction="apca", dims=1:8,
                        reduction.name="adt.umap")
}
# Elbow plot for ADT PCs
print(ElbowPlot(merged_obj, reduction = "apca", ndims = 40))

# 6) Cell-type annotation (after integration)


# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")
# merged_obj <- run_sctype(
#       merged_obj, assay="RNA", scaled=TRUE,
#       known_tissue_type="Immune system",
#       custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",
#       name="sctype_classification"
#     )
  


# Save
saveRDS(merged_obj, save_rds)
cat("Saved:", save_rds, "\n")


# 7) VDJ summaries (only if columns exist)
if ("t_cdr3s_aa" %in% colnames(merged_obj@meta.data)) {
  tcr <- table(merged_obj$t_cdr3s_aa); tcr <- tcr[!is.na(names(tcr)) & nzchar(names(tcr))]
  cat("Top TCR:\n"); print(head(sort(tcr, TRUE), 10))
}
if ("b_cdr3s_aa" %in% colnames(merged_obj@meta.data)) {
  bcr <- table(merged_obj$b_cdr3s_aa); bcr <- bcr[!is.na(names(bcr)) & nzchar(names(bcr))]
  cat("Top BCR:\n"); print(head(sort(bcr, TRUE), 10))
}

# Core plots
print(DimPlot(merged_obj, reduction="rna.umap", group.by="batch"))
print(DimPlot(merged_obj, reduction="adt.umap", group.by="batch"))

if ("sctype_classification" %in% colnames(merged_obj@meta.data)) {
  print(DimPlot(merged_obj, reduction="rna.umap", label=FALSE, repel=TRUE,
                group.by="sctype_classification", label.size=6))
}


if ("sctype_classification" %in% colnames(merged_obj@meta.data)) {
  print(DimPlot(merged_obj, reduction="adt.umap", label=FALSE, repel=TRUE,
                group.by="sctype_classification", label.size=6))
}




VlnPlot(merged_obj, features = "nCount_ADT", group.by = "batch", assay="ADT")
RidgePlot(merged_obj, features = rownames(merged_obj[["ADT"]])[1:8],
          group.by="batch", assay="ADT")



 