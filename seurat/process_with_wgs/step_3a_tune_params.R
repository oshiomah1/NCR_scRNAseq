

suppressPackageStartupMessages({
  pkgs <- c(
    "Seurat", "harmony", "ggplot2", "dplyr",
    "scales", "readxl", "stringr", "tibble",
    "openxlsx", "patchwork"
  )
  print(sapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "MISSING")))
  invisible(lapply(pkgs, require, character.only = TRUE))
  library(future)
})
library(HGNChelper)

plan("sequential")
options(future.globals.maxSize = Inf)
set.seed(411)

## ------------------------------------------------------------------------
## Paths (EDIT THESE IF NEEDED)
## ------------------------------------------------------------------------

log_file   <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/harmony_sctype_tester_log.txt"
output_pdf <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/harmony_sctype_tester_plots.pdf"

# NOTE: you currently override the 5-batch list with only batch5 + batch3.
# Comment out the second 'paths <-' if you want all 5 batches.

paths <- c(
  "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch1_paths_combined_wdj.rds",
  "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch2_paths_combined_wdj.rds",
  "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch3_paths_combined_wdj.rds",
  "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch4_paths_combined_wdj.rds",
  "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/2_merged_lanes_wgs/batch5_paths_combined_wdj.rds"
)

## ------------------------------------------------------------------------
## Logging / PDF
## ------------------------------------------------------------------------

sink(log_file, split = TRUE)
on.exit({ try(sink(), silent = TRUE) }, add = TRUE)

pdf(output_pdf, width = 8, height = 6)
on.exit({ if (dev.cur() > 1) try(dev.off(), silent = TRUE) }, add = TRUE)

cat("=== Harmony + scType tester script ===\n")

## ------------------------------------------------------------------------
## 1) Load & merge
## ------------------------------------------------------------------------

seurat_list <- lapply(paths, readRDS)
merged_obj  <- Reduce(function(x, y) merge(x, y), seurat_list)
rm(seurat_list); gc()

DefaultAssay(merged_obj) <- "RNA"

cat("Cells:", ncol(merged_obj), "\n")
cat("Batches:\n")
print(table(merged_obj$batch))

## ------------------------------------------------------------------------
## 2) RNA preprocessing + PCA
## ------------------------------------------------------------------------

merged_obj <- NormalizeData(merged_obj, normalization.method = "LogNormalize", verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 100, reduction.name = "pca", reduction.key = "PCA_", verbose = FALSE)

cat("Finished PCA.\n")

## PCA diagnostics
cat("Plotting PCA diagnostics (ElbowPlot + DimHeatmaps)...\n")
print(ElbowPlot(merged_obj, ndims = 60))

#print(DimHeatmap(merged_obj, dims = 1:12, cells = 500, balanced = TRUE))
#print(DimHeatmap(merged_obj, dims = 13:24, cells = 500, balanced = TRUE))

## ------------------------------------------------------------------------
## 3) Run Harmony on max PCs
## ------------------------------------------------------------------------

dims_max <- 100   # max PCs to give Harmony
cat("Running Harmony on dims 1:", dims_max, "...\n")

merged_obj <- RunHarmony(
  object = merged_obj,
  group.by.vars = "batch",
  dims.use = 1:dims_max,
  assay.use = "RNA",
  plot_convergence = TRUE
)

cat("Harmony completed.\n")

## Pre vs post batch mixing
p_pca_batch <- DimPlot(merged_obj, reduction = "pca", group.by = "batch") +
  ggtitle("PCA colored by batch")

p_harmony_batch <- DimPlot(merged_obj, reduction = "harmony", group.by = "batch") +
  ggtitle("Harmony colored by batch")

print(p_pca_batch)
print(p_harmony_batch)

## ------------------------------------------------------------------------
## 4) Run scType ONCE (independent of Harmony dims/res)
## ------------------------------------------------------------------------

cat("Running scType annotation...\n")


baseline_dims <- 50
baseline_res  <- 0.6

merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:baseline_dims)
merged_obj <- FindClusters(merged_obj, resolution = baseline_res)
table(merged_obj$seurat_clusters)
# Local paths to scType files
sctype_wrapper_path <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/sctype_wrapper.R"
sctype_db_path      <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/ScTypeDB_short.xlsx"
print(sctype_wrapper_path)

sctype_success <- FALSE

if (file.exists(sctype_wrapper_path)) {
  cat("Sourcing local scType wrapper from:\n  ", sctype_wrapper_path, "\n")
  source(sctype_wrapper_path)
  sctype_success <- TRUE
} else {
  cat("Local scType wrapper not found at:\n  ", sctype_wrapper_path, "\nTrying to fetch from GitHub (may fail on cluster)...\n")
  tryCatch({
    options(timeout = max(600, getOption("timeout")))
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")
    sctype_success <- TRUE
  }, error = function(e) {
    cat("scType download failed, will skip scType-related steps.\nError:\n", conditionMessage(e), "\n")
  })
}

if (sctype_success) {
  merged_obj <- run_sctype(
    merged_obj,
    assay = "RNA",
    scaled = TRUE,  # we already ran ScaleData on RNA
    known_tissue_type = "Immune system",
    custom_marker_file = sctype_db_path,
    name = "sctype_classification"
  )
  
  cat("scType annotation done. Example label counts:\n")
  print(head(table(merged_obj$sctype_classification)))
  
  ## Quick UMAP to visualize scType labels (will be overwritten later)
  merged_obj <- RunUMAP(
    merged_obj,
    reduction      = "harmony",
    dims           = 1:50,
    reduction.name = "rna.umap",
    reduction.key  = "rnaUMAP_",
    verbose        = FALSE
  )
  
  print(
    DimPlot(
      merged_obj,
      reduction = "rna.umap",
      group.by  = "sctype_classification",
      label     = TRUE,
      repel     = TRUE
    ) + ggtitle("UMAP – scType labels (dims 1:50)")
  )
} else {
  cat("Skipping scType annotation and scType-based plots because wrapper could not be loaded.\n")
}

## ------------------------------------------------------------------------
## 5) Grid over dims & resolutions:
##    neighbor graph, clustering, UMAP, plots
## ------------------------------------------------------------------------

dims_grid <- c(40,50, 60)          # candidate Harmony dims
res_grid  <- c(0.6, 0.8, 1.0, 1.2)     # candidate clustering resolutions

cat("Exploring dims_grid =", paste(dims_grid, collapse = ", "),
    "and res_grid =", paste(res_grid, collapse = ", "), "\n")

for (d in dims_grid) {
  cat("\n=== Using Harmony dims 1:", d, "===\n")
  
  ## IMPORTANT: build neighbors on HARMONY, not on UMAP
  merged_obj <- FindNeighbors(
    merged_obj,
    reduction = "harmony",
    dims      = 1:d,
    verbose   = FALSE
  )
  
  ## Compute a UMAP for this 'd'
  merged_obj <- RunUMAP(
    merged_obj,
    reduction      = "harmony",
    dims           = 1:d,
    reduction.name = "rna.umap",
    reduction.key  = "rnaUMAP_",
    verbose        = FALSE
  )
  
  for (res in res_grid) {
    cat("  -> Resolution:", res, "\n")
    
    merged_obj <- FindClusters(merged_obj, resolution = res, verbose = FALSE)
    res_col <- paste0("RNA_snn_res.", res)
    
    # UMAP by clusters
    p_clusters <- DimPlot(
      merged_obj,
      reduction = "rna.umap",
      group.by  = res_col,
      label     = TRUE,
      repel     = TRUE
    ) +
      ggtitle(paste0("UMAP – clusters (dims 1:", d, ", res=", res, ")"))
    
    # UMAP by batch
    p_batch <- DimPlot(
      merged_obj,
      reduction = "rna.umap",
      group.by  = "batch"
    ) +
      ggtitle(paste0("UMAP – batch (dims 1:", d, ", res=", res, ")"))
    
    print(p_clusters)
    
    # Optional: UMAP by scType on the same embedding
    if ("sctype_classification" %in% colnames(merged_obj@meta.data)) {
      p_sctype <- DimPlot(
        merged_obj,
        reduction = "rna.umap",
        group.by  = "sctype_classification",
        label     = TRUE,
        repel     = TRUE
      ) + ggtitle(paste0("UMAP – scType (dims 1:", d, ", res=", res, ")"))
      
      print(p_sctype)
    }
  }
}

## ------------------------------------------------------------------------
## 6) Helper: evaluate how scType labels align with clusters
## ------------------------------------------------------------------------

if (!"sctype_classification" %in% colnames(merged_obj@meta.data)) {
  cat("sctype_classification not found; skipping confusion tables.\n")
} else {
  # Find all clustering resolutions you’ve computed
  res_cols <- grep("^RNA_snn_res\\.", colnames(merged_obj@meta.data), value = TRUE)
  
  cat("Found clustering resolutions for confusion tables:\n")
  print(res_cols)
  
  for (res_col in res_cols) {
    res_val <- sub("RNA_snn_res\\.", "", res_col)
    
    cat("\n==============================\n")
    cat("Resolution:", res_val, " (column:", res_col, ")\n")
    cat("==============================\n")
    
    clust_vec  <- merged_obj[[res_col, drop = TRUE]]
    sctype_vec <- merged_obj$sctype_classification
    
    # Raw counts: clusters x scType labels
    tab <- table(cluster = clust_vec, sctype = sctype_vec)
    cat("\nCluster x scType counts:\n")
    print(tab)
    
    # Row-wise proportions (cluster purity)
    cat("\nRow-wise proportions (per cluster):\n")
    print(round(prop.table(tab, margin = 1), 3))
  }
  
  ## ----------------------------------------------------------------------
  ## 7) Optional: UMAPs for each resolution vs scType (using last UMAP)
  ## ----------------------------------------------------------------------
  
  cat("\nPlotting UMAPs for each resolution vs scType...\n")
  
  for (res_col in res_cols) {
    Idents(merged_obj) <- merged_obj[[res_col, drop = TRUE]]
    
    p_clusters <- DimPlot(
      merged_obj,
      reduction = "rna.umap",
      group.by  = res_col,
      label     = TRUE,
      repel     = TRUE
    ) +
      ggtitle(paste0("Clusters (", res_col, ")"))
    
    p_sctype <- DimPlot(
      merged_obj,
      reduction = "rna.umap",
      group.by  = "sctype_classification",
      label     = TRUE,
      repel     = TRUE
    ) +
      ggtitle("scType labels")
    
    print(p_clusters)
    print(p_sctype)
  }
}

## ------------------------------------------------------------------------
## 8) Optional: clustree across resolutions
## ------------------------------------------------------------------------

if (requireNamespace("clustree", quietly = TRUE)) {
  cat("\nMaking clustree plot for all RNA_snn_res.* columns...\n")
  library(clustree)
  print(clustree::clustree(merged_obj, prefix = "RNA_snn_res."))
} else {
  cat("\nPackage 'clustree' not installed; skipping clustree plot.\n")
}

cat("\n=== Tester script finished. Plots written to:\n", output_pdf, "\n===\n")
