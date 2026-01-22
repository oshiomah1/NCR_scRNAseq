library(Seurat)
library(dplyr)
#library(purrr)
#library(tidyr)
library(ggplot2)
library(tidyverse)
#seur_obj <-readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/raw_merge_all_batches_harm_sc_annotated_all_res6_pca15.rds")
seur_obj <-readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/ScType_multiDB_out_res12_pca15_oscar/raw_merge_all_batches_harm_annotated_all_res12_pca20_noann_Seurat_ScType_6DB_oscar.rds")

 

#read in meta data from suerat object and select relevant columns
meta <- seur_obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%   # add barcodes
  dplyr::select(
    cell,
    NCR.ID,
    validated_TB_status,
    Age.,
    Gender,
    sctype_default, group
  ) %>%
  dplyr::filter(!is.na(NCR.ID))

unique(meta$sctype_default)

meta$celltype_simplified <- dplyr::case_when(
  # CD4
  meta$sctype_default %in% c("Naive CD4+ T cells", "Memory CD4+ T cells") ~ "CD4 T cells",
  
  # CD8
  meta$sctype_default %in% c("Naive CD8+ T cells", "Effector CD8+ T cells") ~ "CD8 T cells",
  
  # NKT-like (often TCR+ cytotoxic; keep separate from NK/CD8 if you're tracking it)
  meta$sctype_default %in% c("CD8+ NKT-like cells") ~ "NKT-like",
  
  # NK (note the double-space in your label)
  meta$sctype_default %in% c("Natural killer  cells") ~ "NK cells",
  
  # gamma delta
  meta$sctype_default %in% c("γδ-T cells") ~ "Gamma delta T cells",
  
  # B
  meta$sctype_default %in% c("Naive B cells", "Plasma B cells") ~ "B cells", # 
  
  # Myeloid: monocytes
  meta$sctype_default %in% c("Classical Monocytes", "Non-classical monocytes") ~ "Monocytes",
  
  # DCs
  meta$sctype_default %in% c("Myeloid Dendritic cells", "Plasmacytoid Dendritic cells") ~ "Dendritic cells",
  
  # Rare types
  meta$sctype_default %in% c("Megakaryocyte") ~ "Megakaryocytes",
  meta$sctype_default %in% c("Mast cells") ~ "Mast cells",
  
  # ISG program
  meta$sctype_default %in% c("ISG expressing immune cells") ~ "ISG-high",
  
  # Unknown
  meta$sctype_default %in% c("Unknown") ~ "Unknown",
  
  TRUE ~ meta$sctype_default
)

# write back to object
#make sure seurat object has this
seur_obj$celltype_simplified <- meta$celltype_simplified[
  match(colnames(seur_obj), meta$cell)
]
unique(meta$celltype_simplified)
table(seur_obj$celltype_simplified, useNA = "ifany")
#join count layers across the batches
seur_obj[["RNA"]] <- JoinLayers(
  object = seur_obj[["RNA"]],
  layers = "counts"
)



###pseudobulk function made more flexible (new)
pseudobulk_rna_by_meta <- function(
    obj,
    meta_col,
    meta_value,
    donor_col = "NCR.ID"
) {
  
  stopifnot(meta_col %in% colnames(obj@meta.data))
  
  cells_use <- WhichCells(
    obj,
    expression = .data[[meta_col]] == meta_value
  )
  
  if (length(cells_use) == 0) {
    stop("No cells found for ", meta_col, " = ", meta_value)
  }
  
  counts <- GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  )[, cells_use, drop = FALSE]
  
  donors <- factor(obj@meta.data[cells_use, donor_col])
  
  donor_mm <- Matrix::sparse.model.matrix(~ donors - 1)
  colnames(donor_mm) <- levels(donors)
  
  counts %*% donor_mm
}

pseudobulk_all <- function(obj, donor_col = "NCR.ID") {
  
  counts <- GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  )
  
  donors <- factor(obj@meta.data[[donor_col]])
  donor_mm <- Matrix::sparse.model.matrix(~ donors - 1)
  colnames(donor_mm) <- levels(donors)
  
  counts %*% donor_mm
}

min_cells <- 500  # 50 is OK as a starting point for global; many people use 100–500 depending on dataset


donor_counts <- seur_obj@meta.data %>%
  dplyr::filter(!is.na(NCR.ID)) %>%
  dplyr::count(NCR.ID, name = "n_cells")

keep_donors <- donor_counts %>%
  filter(n_cells >= min_cells) %>%
  pull(NCR.ID)

obj_filt <- subset(
  seur_obj,
  cells = WhichCells(seur_obj, expression = NCR.ID %in% keep_donors)
)

pb_all_cells <- pseudobulk_all(obj_filt)


donor_meta <- obj_filt@meta.data %>%
  dplyr::filter(!is.na(NCR.ID)) %>%
  dplyr::select(NCR.ID, validated_TB_status, Age., Gender, group) %>%
  dplyr::distinct()


library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)



run_deseq_one_celltype <- function(
    counts,
    celltype,
    donor_meta,
    status_col = "validated_TB_status",
    ref_level = "ctrl",
    case_level = "2weekCase",
    design_formula = ~ Age_c + Gender + group + validated_TB_status,
    min_donors_per_group = 2,
    top_n = 50,
    outdir = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/deseq_res12_15_default"
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- coldata aligned to count columns ----
  coldata <- donor_meta[donor_meta$NCR.ID %in% colnames(counts), ]
  rownames(coldata) <- coldata$NCR.ID
  coldata <- coldata[colnames(counts), , drop = FALSE]
  
  if (!all(rownames(coldata) == colnames(counts))) {
    stop("Donor metadata not aligned for celltype: ", celltype)
  }
  
  # covariates
  coldata$Gender <- factor(coldata$Gender)
  coldata$group  <- factor(coldata$group)
  coldata[[status_col]] <- factor(coldata[[status_col]])
  coldata$Age_c <- scale(coldata$Age., center = TRUE, scale = FALSE)
  
  # set reference level
  coldata[[status_col]] <- relevel(coldata[[status_col]], ref = ref_level)
  
  # check group sizes
  tab <- table(coldata[[status_col]])
  if (!(case_level %in% names(tab)) || !(ref_level %in% names(tab))) {
    message("Skipping ", celltype, ": missing case or ref level.")
    return(NULL)
  }
  if (tab[[case_level]] < min_donors_per_group || tab[[ref_level]] < min_donors_per_group) {
    message("Skipping ", celltype, ": too few donors per group (",
            case_level, "=", tab[[case_level]], ", ",
            ref_level, "=", tab[[ref_level]], ").")
    return(NULL)
  }
  
  # Ensure integer matrix for DESeq2
  counts_mat <- as.matrix(counts)
  storage.mode(counts_mat) <- "numeric"
  counts_mat <- round(counts_mat)
  storage.mode(counts_mat) <- "integer"
  
  # ---- DESeq2 ----
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = coldata,
    design = design_formula
  )
  
  # gene filter
  #dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 10, ]
  dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 70, ]
  dds <- DESeq2::DESeq(dds)
  
  res <- DESeq2::results(dds, contrast = c(status_col, case_level, ref_level))
  res_ordered <- res[order(res$padj), ]
  
  # shrink LFC for the status effect
  coef_name <- grep(paste0("^", status_col), DESeq2::resultsNames(dds), value = TRUE)[1]
  resLFC <- DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  # combine: padj from res, shrunken LFC from resLFC
  plot_df <- as.data.frame(res)
  plot_df$log2FoldChange <- resLFC$log2FoldChange
  
  # ---- Save tables ----
  safe_ct <- gsub("[^A-Za-z0-9_]+", "_", celltype)
  
  write.csv(as.data.frame(res_ordered),
            file = file.path(outdir, paste0("DESeq2_", safe_ct, "_res_ordered.csv")))
  write.csv(as.data.frame(plot_df),
            file = file.path(outdir, paste0("DESeq2_", safe_ct, "_res_with_shrunkLFC.csv")))
  
  # ---- Volcano ----
  png(file.path(outdir, paste0("Volcano_", safe_ct, ".png")),
      width = 1400, height = 900, res = 150)
  print(EnhancedVolcano::EnhancedVolcano(
    plot_df,
    lab = rownames(plot_df),
    x = "log2FoldChange",
    y = "padj",
    title = paste0(celltype, " (", case_level, " vs ", ref_level, ")")
  ))
  dev.off()
  
  # ---- Heatmaps ----
  top_genes <- rownames(head(res_ordered[!is.na(res_ordered$padj), ], top_n))
  
  if (length(top_genes) >= 2) {
    vsd <- DESeq2::vst(dds, blind = FALSE)
    mat <- SummarizedExperiment::assay(vsd)[top_genes, , drop = FALSE]
    mat_z <- t(scale(t(mat)))
    
    ann_col <- as.data.frame(SummarizedExperiment::colData(dds)[, c(status_col, "Gender"), drop = FALSE])
    ann_col[[status_col]] <- factor(ann_col[[status_col]], levels = c(ref_level, case_level))
    ann_col$Gender <- factor(ann_col$Gender)
    
    # clustered columns
    pheatmap::pheatmap(
      mat_z,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      annotation_col = ann_col,
      show_colnames = FALSE,
      fontsize_row = 8,
      main = paste0("Top ", top_n, " DE genes (", celltype, ")"),
      filename = file.path(outdir, paste0("Heatmap_clustered_", safe_ct, ".png"))
    )
    
    # ordered by status
    ord <- order(ann_col[[status_col]])
    mat_z2 <- mat_z[, ord, drop = FALSE]
    ann2 <- ann_col[ord, , drop = FALSE]
    
    pheatmap::pheatmap(
      mat_z2,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_col = ann2,
      show_colnames = FALSE,
      fontsize_row = 8,
      main = paste0("Top ", top_n, " DE genes (", celltype, ") — ordered by status"),
      filename = file.path(outdir, paste0("Heatmap_byStatus_", safe_ct, ".png"))
    )
  } else {
    message("Skipping heatmaps for ", celltype, ": too few valid genes for a heatmap.")
  }
  
  # return objects
  list(dds = dds, res = res, res_ordered = res_ordered, resLFC = resLFC)
}


res_global <- run_deseq_one_celltype(
  counts = pb_all_cells,
  celltype = "Global_PBMC", #label
  donor_meta = donor_meta,
  outdir = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/deseq_global_20pcs_strict10_70"
)

library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(ggplot2)

# 1) Pathways (Hallmark)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

# 2) DESeq2 results -> data.frame
res_df <- as.data.frame(res_global$res)

# 3) Make ranks (uses your function)
ranks <- make_rank(res_df, rank_col = "stat")

# 4) Run fgsea
fg_global <- run_fgsea_one(ranks, pathways = pathways, minSize = 15, maxSize = 500)

# 5) Save
saveRDS(fg_global, "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/GSEA_results/GSEA_global_hallmark.rds")
write.csv(fg_global, "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/GSEA_results/GSEA_global_hallmark.csv",
          row.names = FALSE)


# Choose top pathways overall (by best padj across any cell type)
# top_pathways_global2 <- fg_global %>%
#   filter(!is.na(padj)) %>%
#   group_by(pathway) %>%
#   summarise(best_padj = min(padj), .groups = "drop") %>%
#   arrange(best_padj) %>%
#   slice_head(n = 30) %>%
#   pull(pathway)
# 
# plot_df_global2 <- fg_global %>%
#   filter(pathway %in% top_pathways_global) %>%
#   mutate(
#     neglog10_fdr = -log10(padj + 1e-300),
#     pathway = factor(pathway, levels = rev(top_pathways_global))
#   )


top_df_global <- fg_global %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 30) %>%
  mutate(pathway = factor(pathway, levels = rev(pathway)))

gsea_global2 <- ggplot(top_df, aes(x = pathway, y = NES)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(title = "Top GSEA pathways (Global PBMC pseudobulk)", x = NULL, y = "NES")

gsea_global2

ggsave(
  filename = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/GSEA_results/GSEA_global_hallmarkonly.png",
  plot = gsea_global2,
  width = 12, height = 8, units = "in", dpi = 300
)


library(dplyr)
 

# top_pathways_global <- fg_global %>%
#   filter(!is.na(padj)) %>%
#   arrange(padj) %>%
#   slice_head(n = 30) %>%
#   pull(pathway)
# 
# plot_df_global <- fg_global %>%
#   filter(pathway %in% top_pathways_global) %>%
#   mutate(
#     celltype = "Global_PBMC",
#     neglog10_fdr = -log10(padj + 1e-300),
#     pathway = factor(pathway, levels = rev(top_pathways_global))
#   )
# 
# p <- ggplot(plot_df_global, aes(x = celltype, y = pathway)) +
#   geom_point(aes(size = neglog10_fdr, color = NES)) +
#   theme_bw(base_size = 11) +
#   labs(
#     title = "GSEA (Global PBMC pseudobulk)",
#     x = NULL, y = NULL,
#     size = "-log10(FDR)",
#     color = "NES"
#   )
# 
# p


# GO via MSigDB (C5): BP / MF / CC
# --- Gene Ontology collections from MSigDB (C5) ---
go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  select(gs_name, gene_symbol)
go_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  select(gs_name, gene_symbol)
go_cc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>%
  select(gs_name, gene_symbol)

pathways_go_bp <- split(go_bp$gene_symbol, go_bp$gs_name)
pathways_go_mf <- split(go_mf$gene_symbol, go_mf$gs_name)
pathways_go_cc <- split(go_cc$gene_symbol, go_cc$gs_name)

#   run separately  - i'm using ranks made from deseq
fg_go_bp <- run_fgsea_one(ranks, pathways_go_bp, minSize = 15, maxSize = 500) %>%
  mutate(ontology = "BP")
fg_go_mf <- run_fgsea_one(ranks, pathways_go_mf, minSize = 15, maxSize = 500) %>%
  mutate(ontology = "MF")
fg_go_cc <- run_fgsea_one(ranks, pathways_go_cc, minSize = 15, maxSize = 500) %>%
  mutate(ontology = "CC")

fg_go_all <- bind_rows(fg_go_bp, fg_go_mf, fg_go_cc)

 

top_go <- fg_go_all %>%
  filter(!is.na(padj)) %>%
  group_by(ontology) %>%
  arrange(padj) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  mutate(term = factor(pathway, levels = rev(unique(pathway))))

p_go <- ggplot(top_go, aes(x = term, y = NES)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ ontology, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(title = "GO GSEA (Global PBMC pseudobulk)", x = NULL, y = "NES")

p_go

ggsave(
  filename = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/GSEA_results/global_GO.png",
  plot = p_go,
  width = 22, height = 8, units = "in", dpi = 300
)


