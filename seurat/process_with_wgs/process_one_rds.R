args <- commandArgs(trailingOnly = TRUE)
rds <- args[1]
outdir_base <- args[2]

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(openxlsx); library(patchwork) ;library(openxlsx)
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

tcell_rna_markers <- na.omit(markers[markers$Cell.type == "T-cells", "RNA_name"]) |> as.character()
bcell_rna_markers <- na.omit(markers[markers$Cell.type == "B-cell", "RNA_name"]) |> as.character()
monocytes_rna_markers <- na.omit(markers[markers$Cell.type == "Monocytes", "RNA_name"]) |> as.character()
nkcell_rna_markers <- na.omit(markers[markers$Cell.type == "NK-cells", "RNA_name"]) |> as.character()
dendritic_cell_rna_markers <- na.omit(markers[markers$Cell.type == "Dendritic-Cell", "RNA_name"]) |> as.character()
adt_markers <- na.omit(markers$ADT_name) |> as.character()

# --- ScType setup ---
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/sctrype_noTCR.xlsx"
tissue <- "Immune system"
gs_list <- gene_sets_prepare(db_, tissue)

# --- helpers  ---
safe_ggsave <- function(filename, plot, width=14, height=7, dpi=300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = filename, plot = plot, width = width, height = height, units = "in", dpi = dpi)
}

get_scaled_matrix <- function(seur_obj, assay = "RNA") {
  seurat_package_v5 <- isFALSE("counts" %in% names(attributes(seur_obj[[assay]])))
  if (seurat_package_v5) as.matrix(seur_obj[[assay]]$scale.data) else as.matrix(seur_obj[[assay]]@scale.data)
}

run_sctype_and_label <- function(seur_obj, gs_list) {
  scRNAseqData_scaled <- get_scaled_matrix(seur_obj, assay = "RNA")
  es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE,
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  cl_results <- do.call("rbind", lapply(unique(seur_obj@meta.data$seurat_clusters), function(cl) {
    cells_in_cl <- rownames(seur_obj@meta.data[seur_obj@meta.data$seurat_clusters == cl, , drop = FALSE])
    es.max.cl <- sort(rowSums(es.max[, cells_in_cl, drop = FALSE]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = length(cells_in_cl)), 10)
  }))
  
  sctype_scores <- cl_results %>%
    group_by(cluster) %>%
    slice_max(scores, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  sctype_scores$type[as.numeric(sctype_scores$scores) < (sctype_scores$ncells / 4)] <- "Unknown"
  
  seur_obj@meta.data$sctype_classification_man2 <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, , drop = FALSE]
    seur_obj@meta.data$sctype_classification_man2[seur_obj@meta.data$seurat_clusters == j] <- as.character(cl_type$type[1])
  }
  
  list(seur_obj = seur_obj, cl_results = cl_results)
}

make_barplot_inset <- function(seur_obj, prefix, outdir_bar,
                               celltype_col = "sctype_classification_man",
                               zoom_max = 0.8) {
  dir.create(outdir_bar, recursive = TRUE, showWarnings = FALSE)
  
  # Make sure required columns exist
  needed <- c("NCR.ID", "validated_TB_status", celltype_col)
  missing <- setdiff(needed, colnames(seur_obj@meta.data))
  if (length(missing) > 0) {
    message("Skipping barplot for ", prefix, " (missing: ", paste(missing, collapse=", "), ")")
    return(invisible(NULL))
  }
  
  # Aggregate cell counts per donor
  celltype_props <- seur_obj@meta.data %>%
    dplyr::filter(!is.na(NCR.ID)) %>%
    dplyr::count(NCR.ID, validated_TB_status, .data[[celltype_col]], name = "n_cells") %>%
    dplyr::group_by(NCR.ID) %>%
    dplyr::mutate(
      total_cells = sum(n_cells),
      prop = n_cells / total_cells
    ) %>%
    dplyr::ungroup()
  
  # donor-level metadata (Age, Gender) if present
  donor_cols <- intersect(c("NCR.ID", "Age.", "Gender"), colnames(seur_obj@meta.data))
  donor_meta <- seur_obj@meta.data %>%
    dplyr::filter(!is.na(NCR.ID)) %>%
    dplyr::select(dplyr::all_of(donor_cols)) %>%
    dplyr::distinct()
  
  celltype_props <- celltype_props %>%
    dplyr::left_join(donor_meta, by = "NCR.ID")
  
  # Join antibiotic metadata (optional)
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
  
  # Summary per TB status x cell type
  summary_df <- celltype_props %>%
    dplyr::group_by(validated_TB_status, .data[[celltype_col]]) %>%
    dplyr::summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(percent = mean_prop * 100)
  
  # Order cell types by overall mean
  celltype_order <- summary_df %>%
    dplyr::group_by(.data[[celltype_col]]) %>%
    dplyr::summarise(overall_mean = mean(percent), .groups = "drop") %>%
    dplyr::arrange(overall_mean) %>%
    dplyr::pull(.data[[celltype_col]])
  
  summary_df[[celltype_col]] <- factor(summary_df[[celltype_col]], levels = celltype_order)
  
  # Main plot
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
  
  # Inset for bottom 3 cell types
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

res <- run_sctype_and_label(seur_obj, gs_list)
seur_obj <- res$seur_obj
cl_results <- res$cl_results


# output dirs
out_umap <- file.path(outdir_base, "umap")
out_dot  <- file.path(outdir_base, "dotplots")
out_top4 <- file.path(outdir_base, "top4_per_cluster")
out_bar <- file.path(outdir_base, "cell_type_proportions")

dir.create(out_umap, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dot,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_top4, recursive = TRUE, showWarnings = FALSE)
dir.create(out_bar, recursive = TRUE, showWarnings = FALSE)

# UMAP
p_umap <- DimPlot(seur_obj, reduction="rna.umap", label=TRUE, repel=TRUE, group.by="sctype_classification_man2")
safe_ggsave(file.path(out_umap, paste0(prefix, "_umap_sctype_man2.png")), p_umap, width=14, height=7)


# Dotplots
dp_t <- DotPlot(seur_obj, features=tcell_rna_markers, group.by="sctype_classification_man2", assay="RNA") + RotatedAxis() + ggtitle("T cells (RNA)")
safe_ggsave(file.path(out_dot, paste0(prefix, "_dotplot_Tcell_RNA.png")), dp_t, width=18, height=7)

dp_nk <- DotPlot(seur_obj, features=nkcell_rna_markers, group.by="sctype_classification_man2", assay="RNA") + RotatedAxis() + ggtitle("NK (RNA)")
safe_ggsave(file.path(out_dot, paste0(prefix, "_dotplot_NK_RNA.png")), dp_nk, width=18, height=7)

dp_b <- DotPlot(seur_obj, features=bcell_rna_markers, group.by="sctype_classification_man2", assay="RNA") + RotatedAxis() + ggtitle("B cells (RNA)")
safe_ggsave(file.path(out_dot, paste0(prefix, "_dotplot_Bcell_RNA.png")), dp_b, width=18, height=7)

dp_m <- DotPlot(seur_obj, features=monocytes_rna_markers, group.by="sctype_classification_man2", assay="RNA") + RotatedAxis() + ggtitle("Monocytes (RNA)")
safe_ggsave(file.path(out_dot, paste0(prefix, "_dotplot_Monocytes_RNA.png")), dp_m, width=18, height=7)

dp_dc <- DotPlot(seur_obj, features=dendritic_cell_rna_markers, group.by="sctype_classification_man2", assay="RNA") + RotatedAxis() + ggtitle("DC (RNA)")
safe_ggsave(file.path(out_dot, paste0(prefix, "_dotplot_DC_RNA.png")), dp_dc, width=18, height=7)

if ("ADT" %in% Assays(seur_obj)) {
  dp_adt <- DotPlot(seur_obj, features=adt_markers, group.by="sctype_classification_man2", assay="ADT") + RotatedAxis() + ggtitle("ADT")
  safe_ggsave(file.path(out_dot, paste0(prefix, "_dotplot_ADT.png")), dp_adt, width=18, height=7)
}

# Top4 PDF
top4 <- cl_results %>%
  group_by(cluster) %>%
  slice_max(scores, n=4, with_ties=FALSE) %>%
  ungroup()

pdf(file.path(out_top4, paste0(prefix, "_ScType_top4_per_cluster.pdf")), width=7, height=5)
for (cl in sort(unique(top4$cluster))) {
  df_cl <- dplyr::filter(top4, cluster == cl)
  p <- ggplot(df_cl, aes(x=reorder(type, scores), y=scores, fill=type)) +
    geom_col(show.legend=FALSE) + coord_flip() +
    theme_minimal(base_size=12) +
    labs(title=paste("Cluster", cl, "- Top 4 ScType"), x=NULL, y="ScType score")
  print(p)
}
dev.off()


#make barplots tocpmapre ct props

make_barplot_inset(seur_obj, prefix, out_bar, celltype_col = "sctype_classification_man2")
message("DONE: ", prefix)
