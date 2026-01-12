args <- commandArgs(trailingOnly = TRUE)
rds <- args[1]
outdir_base <- args[2]

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(patchwork)
})

library(HGNChelper)
library(stringr)
library(readxl)
set.seed(411)

# --- markers (loaded per job) ---
markers <- read.csv(
  "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/Markers.csv",
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

tcell_rna_markers       <- na.omit(markers[markers$Cell.type == "T-cells",        "RNA_name"]) |> as.character()
bcell_rna_markers       <- na.omit(markers[markers$Cell.type == "B-cell",         "RNA_name"]) |> as.character()
monocytes_rna_markers   <- na.omit(markers[markers$Cell.type == "Monocytes",      "RNA_name"]) |> as.character()
nkcell_rna_markers      <- na.omit(markers[markers$Cell.type == "NK-cells",       "RNA_name"]) |> as.character()
dendritic_cell_rna_markers <- na.omit(markers[markers$Cell.type == "Dendritic-Cell","RNA_name"]) |> as.character()
adt_markers             <- na.omit(markers$ADT_name) |> as.character()

# --- ScType setup ---
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

tissue <- "Immune system"

# --- helpers ---
safe_ggsave <- function(filename, plot, width=14, height=7, dpi=300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = filename, plot = plot, width = width, height = height, units = "in", dpi = dpi)
}

get_scaled_matrix <- function(seur_obj, assay = "RNA") {
  seurat_package_v5 <- isFALSE("counts" %in% names(attributes(seur_obj[[assay]])))
  if (seurat_package_v5) as.matrix(seur_obj[[assay]]$scale.data) else as.matrix(seur_obj[[assay]]@scale.data)
}

run_sctype_and_label <- function(seur_obj, gs_list, label_col) {
  scRNAseqData_scaled <- get_scaled_matrix(seur_obj, assay = "RNA")
  
  es.max <- sctype_score(
    scRNAseqData = scRNAseqData_scaled,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  
  cl_results <- do.call("rbind", lapply(unique(seur_obj@meta.data$seurat_clusters), function(cl) {
    cells_in_cl <- rownames(seur_obj@meta.data[seur_obj@meta.data$seurat_clusters == cl, , drop = FALSE])
    es.max.cl <- sort(rowSums(es.max[, cells_in_cl, drop = FALSE]), decreasing = TRUE)
    head(data.frame(
      cluster = cl,
      type    = names(es.max.cl),
      scores  = es.max.cl,
      ncells  = length(cells_in_cl)
    ), 10)
  }))
  
  sctype_scores <- cl_results %>%
    group_by(cluster) %>%
    slice_max(scores, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  sctype_scores$type[as.numeric(sctype_scores$scores) < (sctype_scores$ncells / 4)] <- "Unknown"
  
  # write into a DB-specific metadata column
  seur_obj@meta.data[[label_col]] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, , drop = FALSE]
    seur_obj@meta.data[[label_col]][seur_obj@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
  }
  
  list(seur_obj = seur_obj, cl_results = cl_results, sctype_scores = sctype_scores)
}

make_barplot_inset <- function(seur_obj, prefix, outdir_bar,
                               celltype_col,
                               zoom_max = 0.8) {
  dir.create(outdir_bar, recursive = TRUE, showWarnings = FALSE)
  
  needed <- c("NCR.ID", "validated_TB_status", celltype_col)
  missing <- setdiff(needed, colnames(seur_obj@meta.data))
  if (length(missing) > 0) {
    message("Skipping barplot for ", prefix, " (missing: ", paste(missing, collapse=", "), ")")
    return(invisible(NULL))
  }
  
  celltype_props <- seur_obj@meta.data %>%
    dplyr::filter(!is.na(NCR.ID)) %>%
    dplyr::count(NCR.ID, validated_TB_status, .data[[celltype_col]], name = "n_cells") %>%
    dplyr::group_by(NCR.ID) %>%
    dplyr::mutate(
      total_cells = sum(n_cells),
      prop = n_cells / total_cells
    ) %>%
    dplyr::ungroup()
  
  donor_cols <- intersect(c("NCR.ID", "Age.", "Gender"), colnames(seur_obj@meta.data))
  donor_meta <- seur_obj@meta.data %>%
    dplyr::filter(!is.na(NCR.ID)) %>%
    dplyr::select(dplyr::all_of(donor_cols)) %>%
    dplyr::distinct()
  
  celltype_props <- celltype_props %>%
    dplyr::left_join(donor_meta, by = "NCR.ID")
  
  extra_meta_path <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/extra_metadata.tsv"
  if (file.exists(extra_meta_path)) {
    abx <- read.csv(extra_meta_path) %>%
      dplyr::rename(days_antibiotics = lab_date_2_trtment_strt) %>%
      dplyr::mutate(
        validated_case_status = factor(validated_case_status),
        days_antibiotics = as.numeric(days_antibiotics)
      )
    celltype_props <- celltype_props %>% dplyr::left_join(abx, by = "NCR.ID")
  } else {
    message("extra_metadata.tsv not found; continuing without ABX join for ", prefix)
  }
  
  summary_df <- celltype_props %>%
    dplyr::group_by(validated_TB_status, .data[[celltype_col]]) %>%
    dplyr::summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(percent = mean_prop * 100)
  
  celltype_order <- summary_df %>%
    dplyr::group_by(.data[[celltype_col]]) %>%
    dplyr::summarise(overall_mean = mean(percent), .groups = "drop") %>%
    dplyr::arrange(overall_mean) %>%
    dplyr::pull(.data[[celltype_col]])
  
  summary_df[[celltype_col]] <- factor(summary_df[[celltype_col]], levels = celltype_order)
  
  p_main <- ggplot(
    summary_df,
    aes(
      x = percent,
      y = .data[[celltype_col]],
      fill = validated_TB_status
    )
  ) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(
      values = c("2weekCase" = "#FFC20A", "ctrl" = "#0C7BDC"),
      name = "TB status"
    ) +
    scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    ) +
    labs(x = "Cell type proportion (%)", y = NULL)
  
  bottom_types <- head(levels(summary_df[[celltype_col]]), 3)
  summary_df_small <- summary_df %>% dplyr::filter(.data[[celltype_col]] %in% bottom_types)
  
  p_inset <- ggplot(
    summary_df_small,
    aes(
      x = percent,
      y = .data[[celltype_col]],
      fill = validated_TB_status
    )
  ) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(
      values = c("2weekCase" = "#FFC20A", "ctrl" = "#0C7BDC"),
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(0, zoom_max),
      breaks = seq(0, zoom_max, by = zoom_max/4),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.text = element_text(size = 9),
      axis.title.x = element_text(size = 9)
    ) +
    labs(x = "Zoomed (%)", y = NULL)
  
  barplots_celltypes_inset <- p_main + patchwork::inset_element(
    p_inset,
    left = 0.45, bottom = 0.08, right = 0.95, top = 0.50
  )
  
  ggsave(
    filename = file.path(outdir_bar, paste0(prefix, "_barplots_celltypes_inset.png")),
    plot = barplots_celltypes_inset,
    width = 18, height = 6, units = "in", dpi = 300
  )
}

# --- run ---
prefix <- sub("\\.rds$", "", basename(rds))
seur_obj <- readRDS(rds)

# output base dirs (we'll create per-DB subfolders underneath each)
out_umap <- file.path(outdir_base, "umap")
out_dot  <- file.path(outdir_base, "dotplots")
out_top4 <- file.path(outdir_base, "top4_per_cluster")
out_bar  <- file.path(outdir_base, "cell_type_proportions")

dir.create(out_umap, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dot,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_top4, recursive = TRUE, showWarnings = FALSE)
dir.create(out_bar,  recursive = TRUE, showWarnings = FALSE)

# --- 7 ScType databases (label -> path) --- #add a no tcr only
dbs <- list(
  baseline                     = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/ScTypeDB_immune_only.xlsx",
  nktlikenegedit        = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/notcr/no_isg_nktlikenegedit.xlsx",
  no_isg                       = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/notcr/no_isg.xlsx",
  no_mega                      = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/notcr/no_mega.xlsx",
  no_isg_no_mega               = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/notcr/no_isg_no_mega.xlsx",
  no_isg_nktlikenegedit        = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/notcr/no_isg_nktlikenegedit.xlsx",
  no_isg_nktlikepozedit        = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/notcr/no_isg_nktlikepozedit.xlsx",
  no_isg_nktlikenegedit_nomega = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/notcr/no_isg_nktlikenegedit_nomega.xlsx"
)

# store cluster-level results in @misc so the final Seurat has everything
if (is.null(seur_obj@misc$sctype_multiDB)) {
  seur_obj@misc$sctype_multiDB <- list()
}
seur_obj@misc$sctype_multiDB$databases <- dbs
seur_obj@misc$sctype_multiDB$cl_results <- list()
seur_obj@misc$sctype_multiDB$top1_per_cluster <- list()

for (db_label in names(dbs)) {
  db_path <- dbs[[db_label]]
  if (!file.exists(db_path)) {
    message("SKIP (missing DB): ", db_label, " -> ", db_path)
    next
  }
  
  message("Running ScType with DB: ", db_label, " (", db_path, ")")
  
  # per-DB subdirs
  out_umap_db <- file.path(out_umap, db_label)
  out_dot_db  <- file.path(out_dot,  db_label)
  out_top4_db <- file.path(out_top4, db_label)
  out_bar_db  <- file.path(out_bar,  db_label)
  
  dir.create(out_umap_db, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dot_db,  recursive = TRUE, showWarnings = FALSE)
  dir.create(out_top4_db, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_bar_db,  recursive = TRUE, showWarnings = FALSE)
  
  # prepare gene sets for this DB
  gs_list <- gene_sets_prepare(db_path, tissue)
  
  # label column name for this DB run (goes into one Seurat object)
  label_col <- paste0("sctype_", db_label)
  
  res <- run_sctype_and_label(seur_obj, gs_list, label_col = label_col)
  seur_obj <- res$seur_obj
  cl_results <- res$cl_results
  sctype_scores <- res$sctype_scores
  
  # save cluster-level outputs into misc
  seur_obj@misc$sctype_multiDB$cl_results[[db_label]] <- cl_results
  seur_obj@misc$sctype_multiDB$top1_per_cluster[[db_label]] <- sctype_scores
  
  # ---- UMAP ----
  p_umap <- DimPlot(
    seur_obj, reduction = "rna.umap",
    label = TRUE, repel = TRUE,
    group.by = label_col
  ) + ggtitle(paste0("ScType — DB: ", db_label))
  safe_ggsave(
    file.path(out_umap_db, paste0(prefix, "__", db_label, "_umap.png")),
    p_umap, width = 14, height = 7
  )
  
  # ---- DotPlots ----
  dp_t <- DotPlot(seur_obj, features = tcell_rna_markers, group.by = label_col, assay = "RNA") +
    RotatedAxis() + ggtitle(paste0("T cells (RNA) — DB: ", db_label))
  safe_ggsave(file.path(out_dot_db, paste0(prefix, "__", db_label, "_dotplot_Tcell_RNA.png")), dp_t, width = 18, height = 7)
  
  dp_nk <- DotPlot(seur_obj, features = nkcell_rna_markers, group.by = label_col, assay = "RNA") +
    RotatedAxis() + ggtitle(paste0("NK (RNA) — DB: ", db_label))
  safe_ggsave(file.path(out_dot_db, paste0(prefix, "__", db_label, "_dotplot_NK_RNA.png")), dp_nk, width = 18, height = 7)
  
  dp_b <- DotPlot(seur_obj, features = bcell_rna_markers, group.by = label_col, assay = "RNA") +
    RotatedAxis() + ggtitle(paste0("B cells (RNA) — DB: ", db_label))
  safe_ggsave(file.path(out_dot_db, paste0(prefix, "__", db_label, "_dotplot_Bcell_RNA.png")), dp_b, width = 18, height = 7)
  
  dp_m <- DotPlot(seur_obj, features = monocytes_rna_markers, group.by = label_col, assay = "RNA") +
    RotatedAxis() + ggtitle(paste0("Monocytes (RNA) — DB: ", db_label))
  safe_ggsave(file.path(out_dot_db, paste0(prefix, "__", db_label, "_dotplot_Monocytes_RNA.png")), dp_m, width = 18, height = 7)
  
  dp_dc <- DotPlot(seur_obj, features = dendritic_cell_rna_markers, group.by = label_col, assay = "RNA") +
    RotatedAxis() + ggtitle(paste0("DC (RNA) — DB: ", db_label))
  safe_ggsave(file.path(out_dot_db, paste0(prefix, "__", db_label, "_dotplot_DC_RNA.png")), dp_dc, width = 18, height = 7)
  
  if ("ADT" %in% Assays(seur_obj)) {
    dp_adt <- DotPlot(seur_obj, features = adt_markers, group.by = label_col, assay = "ADT") +
      RotatedAxis() + ggtitle(paste0("ADT — DB: ", db_label))
    safe_ggsave(file.path(out_dot_db, paste0(prefix, "__", db_label, "_dotplot_ADT.png")), dp_adt, width = 18, height = 7)
  }
  
  # ---- Top4 PDF ----
  top4 <- cl_results %>%
    group_by(cluster) %>%
    slice_max(scores, n = 4, with_ties = FALSE) %>%
    ungroup()
  
  pdf(file.path(out_top4_db, paste0(prefix, "__", db_label, "_ScType_top4_per_cluster.pdf")), width = 7, height = 5)
  for (cl in sort(unique(top4$cluster))) {
    df_cl <- dplyr::filter(top4, cluster == cl)
    p <- ggplot(df_cl, aes(x = reorder(type, scores), y = scores, fill = type)) +
      geom_col(show.legend = FALSE) +
      coord_flip() +
      theme_minimal(base_size = 12) +
      labs(
        title = paste0("Cluster ", cl, " — Top 4 ScType (DB: ", db_label, ")"),
        x = NULL, y = "ScType score"
      )
    print(p)
  }
  dev.off()
  
  # ---- Barplot (cell type proportions) ----
  make_barplot_inset(
    seur_obj,
    prefix = paste0(prefix, "__", db_label),
    outdir_bar = out_bar_db,
    celltype_col = label_col
  )
  
  message("DONE DB: ", db_label)
}

# ---- save ONE final Seurat with all 7 label columns ----
final_rds <- file.path(outdir_base, paste0(prefix, "_Seurat_ScType_7DB.rds"))
saveRDS(seur_obj, final_rds)

message("ALL DONE. Saved combined Seurat: ", final_rds)
