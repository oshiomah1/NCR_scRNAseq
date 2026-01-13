#library(Seurat)
library(dplyr)
#library(purrr)
library(tidyr)
library(ggplot2)
library(tidyverse)
seur_obj <-readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/raw_merge_all_batches_harm_sc_annotated_all_res6_pca15.rds")


#read in meta data from suerat object and select relevant columns
meta <- seur_obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%   # add barcodes
  dplyr::select(
    cell,
    NCR.ID,
    validated_TB_status,
    Age.,
    Gender,
    sctype_classification_man2, group
  ) %>%
  dplyr::filter(!is.na(NCR.ID))

# aggregate sub cell types into broader cell types
# meta$celltype_simplified <- dplyr::case_when(
#   meta$sctype_classification_man %in% c("Naive CD4+ T cells", "Memory CD4+ T cells") ~ "CD4+ T cells",
#   meta$sctype_classification_man %in% c("Naive CD8+ T cells", "Effector CD8+ T cells", "CD8+ NKT-like cells") ~ "CD8+ T cells",
#   meta$sctype_classification_man %in% c("Natural killer  cells") ~ "NK cells",
#   meta$sctype_classification_man %in% c("Naive B cells", "Plasma B cells") ~ "B cells",
#   meta$sctype_classification_man %in% c("Classical Monocytes", "Non-classical monocytes") ~ "Monocytes",
#   meta$sctype_classification_man %in% c("Plasmacytoid Dendritic cells") ~ "Dendritic cells",
#   TRUE ~ meta$sctype_classification_man
# )



meta$celltype_simplified <- dplyr::case_when(
  # CD4
  meta$sctype_classification_man2 %in% c("Naive CD4+ T cells", "Memory CD4+ T cells") ~ "CD4 T cells",
  
  # CD8
  meta$sctype_classification_man2 %in% c("Naive CD8+ T cells", "Effector CD8+ T cells") ~ "CD8 T cells",
  
  # NKT-like (often TCR+ cytotoxic; keep separate from NK/CD8 if you're tracking it)
  meta$sctype_classification_man2 %in% c("CD8+ NKT-like cells") ~ "NKT-like",
  
  # NK (note the double-space in your label)
  meta$sctype_classification_man2 %in% c("Natural killer  cells") ~ "NK cells",
  
  # gamma delta
  meta$sctype_classification_man2 %in% c("Gamma Delta-T cells") ~ "Gamma delta T cells",
  
  # B
  meta$sctype_classification_man2 %in% c("Naive B cells", "Plasma B cells") ~ "B cells", #
  
  # Myeloid: monocytes
  meta$sctype_classification_man2 %in% c("Classical Monocytes", "Non-classical monocytes") ~ "Monocytes",
  
  # DCs
  meta$sctype_classification_man2 %in% c("Myeloid Dendritic cells", "Plasmacytoid Dendritic cells") ~ "Dendritic cells",
  
  # Rare types
  meta$sctype_classification_man2 %in% c("Megakaryocyte") ~ "Megakaryocytes",
  meta$sctype_classification_man2 %in% c("Mast cells") ~ "Mast cells",
  
  # ISG program
  meta$sctype_classification_man2 %in% c("ISG expressing immune cells") ~ "ISG-high",
  
  # Unknown
  meta$sctype_classification_man2 %in% c("Unknown") ~ "Unknown",
  
  TRUE ~ meta$sctype_classification_man2
)

# write back to object
 #make sure seurat object has this
seur_obj$celltype_simplified <- meta$celltype_simplified[
  match(colnames(seur_obj), meta$cell)
]

table(seur_obj$celltype_simplified, useNA = "ifany")
#joimnm count layers across the batches
seur_obj[["RNA"]] <- JoinLayers(
  object = seur_obj[["RNA"]],
  layers = "counts"
)

# pseudobulk function RNA
pseudobulk_rna <- function(obj, celltype, donor_col = "NCR.ID") {
  
  stopifnot("celltype_simplified" %in% colnames(obj@meta.data))
  
  cells_use <- WhichCells(
    obj,
    expression = celltype_simplified == celltype
  )
  
  if (length(cells_use) == 0) {
    stop("No cells found for cell type: ", celltype)
  }
  
  counts <- GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  )[, cells_use, drop = FALSE]
  
  donors <- obj@meta.data[cells_use, donor_col]
  donors <- factor(donors)
  
  #  sparse donor design matrix
  donor_mm <- Matrix::sparse.model.matrix(~ donors - 1)
  colnames(donor_mm) <- levels(donors)
  
  #  pseudobulk aggregation
  pb_counts <- counts %*% donor_mm
  
  return(pb_counts)
}


# EXample usage 
rna_pb_DC <- pseudobulk_rna(seur_obj, "Dendritic cells")

rna_pb_NK <- pseudobulk_rna(seur_obj, "NK cells")


# ADT pseudobullk function
pseudobulk_adt <- function(obj, celltype, donor_col = "NCR.ID") {
  
  stopifnot("celltype_simplified" %in% colnames(obj@meta.data))
  
  cells_use <- WhichCells(
    obj,
    expression = celltype_simplified == celltype
  )
  
  if (length(cells_use) == 0) {
    stop("No cells found for cell type: ", celltype)
  }
  
  counts <- GetAssayData(
    obj,
    assay = "ADT",
    layer = "counts"
  )[, cells_use, drop = FALSE]
  
  donors <- obj@meta.data[cells_use, donor_col]
  donors <- factor(donors)
  
  donor_mm <- Matrix::sparse.model.matrix(~ donors - 1)
  colnames(donor_mm) <- levels(donors)
  
  pb_counts <- counts %*% donor_mm
  
  return(pb_counts)
}

# EXample usage 
adt_pb_NK <- pseudobulk_adt(seur_obj, "NK cells")

#create donor meta

donor_meta <- seur_obj@meta.data %>%
  filter(!is.na(NCR.ID)) %>%
  select(
    NCR.ID,
    validated_TB_status,
    Age.,
    Gender,
    group
  ) %>%
  distinct()


# Replace 'cell_type' with the actual metadata column name in your Seurat object
unique_cell_types <- unique(meta$celltype_simplified)

#now automate the functions across all cell types

# RNA pseudobulk
pseudobulk_rna_results <- lapply(unique_cell_types, function(ct) {
  pseudobulk_rna(seur_obj, ct)
})
names(pseudobulk_rna_results) <- unique_cell_types

# ADT pseudobulk
pseudobulk_adt_results <- lapply(unique_cell_types, function(ct) {
  pseudobulk_adt(seur_obj, ct)
})
names(pseudobulk_adt_results) <- unique_cell_types


saveRDS(pseudobulk_rna_results, file = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/pseudobulk/pseudobulk_rna_results.rds")
saveRDS(pseudobulk_adt_results, file = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/pseudobulk/pseudobulk_adt_results.rds")

#####################################################################
###############READ IN pseudobulk reusults#############################
#####################################################################
#pseudobulk_rna_results<- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/pseudobulk/pseudobulk_rna_results.rds")
#pseudobulk_adt_results <- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/pseudobulk/pseudobulk_adt_results.rds")


#DESeq2 example (RNA, one cell type)
counts <- rna_pb_NK
counts <-adt_pb_NK
library(DESeq2)

# Subset metadata to only donors in counts, in the correct order
coldata <- donor_meta[donor_meta$NCR.ID %in% colnames(counts), ]
rownames(coldata) <- coldata$NCR.ID

# reorder to match count columns
coldata <- coldata[colnames(counts), , drop = FALSE]

# check alignment
all(rownames(coldata) == colnames(counts))  # should be TRUE

# convert covariates to factors
coldata$Gender <- factor(coldata$Gender)
coldata$group <- factor(coldata$group)
coldata$validated_TB_status <- factor(coldata$validated_TB_status)
coldata$Age_c <- scale(coldata$Age., center = TRUE, scale = FALSE) # just subtract mean

#Set ctrl as the reference BEFORE building dds
coldata$validated_TB_status <- relevel(coldata$validated_TB_status, ref = "ctrl")

# create DESeq dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ Age_c + Gender + group + validated_TB_status
)

#remove low expressed genes
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("validated_TB_status", "2weekCase", "ctrl"))

sig_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
res_ordered <- res[order(res$padj), ]
head(res_ordered, 20)

resLFC <- lfcShrink(dds, coef="validated_TB_status_2weekCase_vs_ctrl", type="apeglm")

library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'padj')
library(pheatmap)
#top_genes <- rownames(head(res_ordered, 50))
#pheatmap(assay(vst(dds))[top_genes, ], cluster_rows = TRUE, cluster_cols = TRUE)



# 1) pick top genes (use res, not resLFC, for significance)
top_genes <- rownames(head(res_ordered, 50))

# 2) VST transform (blind=FALSE recommended for DE context)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top_genes, , drop = FALSE]

# 3) Z-score each gene across donors (so color reflects relative up/down)
mat_z <- t(scale(t(mat)))  # row-wise scaling

# 4) Build column annotations from colData(dds)
#ann_col <- as.data.frame(colData(dds)[, c("validated_TB_status","Gender","group","Age_c")])
ann_col <- as.data.frame(colData(dds)[, c("validated_TB_status","Gender")])

ann_col$validated_TB_status <- factor(ann_col$validated_TB_status, levels = c("ctrl","2weekCase"))
ann_col$Gender <- factor(ann_col$Gender)
#ann_col$group <- factor(ann_col$group)

# make Age a numeric column rather than centered scale object
#ann_col$Age_c <- as.numeric(ann_col$Age_c)

pheatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = ann_col,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Top 50 DE genes (VST, row Z-score)"
)

 
ord <- order(ann_col$validated_TB_status)
mat_z2 <- mat_z[, ord]
ann2 <- ann_col[ord, , drop = FALSE]

pheatmap(
  mat_z2,
  cluster_rows = TRUE,
  cluster_cols = FALSE,          # <- makes case/control pattern obvious
  annotation_col = ann2,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Top 50 DE genes (NK) — ordered by TB status"
)





##################
###super functiom###
##################
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
    outdir = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/deseq_res1"
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
  coldata$group <- factor(coldata$group)
  coldata[[status_col]] <- factor(coldata[[status_col]])
  coldata$Age_c <- scale(coldata$Age., center = TRUE, scale = FALSE)
  
  # set reference level
  coldata[[status_col]] <- relevel(coldata[[status_col]], ref = ref_level)
  
  # check group sizes
  tab <- table(coldata[[status_col]])
  if (!(case_level %in% names(tab)) || !(ref_level %in% names(tab))) {
    message("Skipping ", celltype, ": missing case or ref level in this cell type.")
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
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = coldata,
    design = design_formula
  )
  
  dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c(status_col, case_level, ref_level))
  res_ordered <- res[order(res$padj), ]
  
  # shrink LFC for the status effect
  coef_name <- grep(paste0("^", status_col), resultsNames(dds), value = TRUE)[1]
  resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  # combine: padj from res, shrunken LFC from resLFC
  plot_df <- as.data.frame(res)
  plot_df$log2FoldChange <- resLFC$log2FoldChange
  
  # ---- Save tables ----
  safe_ct <- gsub("[^A-Za-z0-9_]+", "_", celltype)
  
  write.csv(
    as.data.frame(res_ordered),
    file = file.path(outdir, paste0("DESeq2_", safe_ct, "_res_ordered.csv"))
  )
  write.csv(
    as.data.frame(plot_df),
    file = file.path(outdir, paste0("DESeq2_", safe_ct, "_res_with_shrunkLFC.csv"))
  )
  
  # ---- Volcano ----
  png(file.path(outdir, paste0("Volcano_", safe_ct, ".png")), width = 1400, height = 900, res = 150)
  print(
    EnhancedVolcano(plot_df,
                    lab = rownames(plot_df),
                    x = "log2FoldChange",
                    y = "padj",
                    title = paste0(celltype, " (", case_level, " vs ", ref_level, ")")
    )
  )
  dev.off()
  
  # ---- Heatmaps ----
  # top genes by padj
  top_genes <- rownames(head(res_ordered[!is.na(res_ordered$padj), ], top_n))
  
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)[top_genes, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  
  ann_col <- as.data.frame(colData(dds)[, c(status_col, "Gender"), drop = FALSE])
  ann_col[[status_col]] <- factor(ann_col[[status_col]], levels = c(ref_level, case_level))
  ann_col$Gender <- factor(ann_col$Gender)
  
  # clustered columns
  png(file.path(outdir, paste0("Heatmap_clustered_", safe_ct, ".png")), width = 1800, height = 900, res = 150)
  pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = ann_col,
    show_colnames = FALSE,
    fontsize_row = 8,
    main = paste0("Top ", top_n, " DE genes (", celltype, ")")
  )
  dev.off()
  
  # ordered by TB status
  ord <- order(ann_col[[status_col]])
  mat_z2 <- mat_z[, ord, drop = FALSE]
  ann2 <- ann_col[ord, , drop = FALSE]
  
  png(file.path(outdir, paste0("Heatmap_byStatus_", safe_ct, ".png")), width = 1800, height = 900, res = 150)
  pheatmap(
    mat_z2,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = ann2,
    show_colnames = FALSE,
    fontsize_row = 8,
    main = paste0("Top ", top_n, " DE genes (", celltype, ") — ordered by status")
  )
  dev.off()
  
  # return objects if you want to inspect in-session
  list(dds = dds, res = res, res_ordered = res_ordered, resLFC = resLFC)
}

# all_results <- lapply(names(pseudobulk_rna_results), function(ct) {
#   run_deseq_one_celltype(
#     counts = pseudobulk_rna_results[[ct]],
#     celltype = ct,
#     donor_meta = donor_meta,
#     outdir = "DESeq2_pseudobulk_by_celltype",
#     top_n = 50,
#     min_donors_per_group = 2
#   )
# })
library(future)
library(future.apply)

# choose cores (don’t take the whole node)
n_workers <- 8

plan(multisession, workers = n_workers)   # or plan(multicore) on Linux if allowed

ct_vector = c("NKT-like","CD8 T cells","NK cells","CD4 T cells", "ISG-high" )               
            
all_results <- future_lapply(names(pseudobulk_rna_results), function(ct) {
  run_deseq_one_celltype(
    counts = pseudobulk_rna_results[[ct]],
    celltype = ct,
    donor_meta = donor_meta,
    outdir = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/DESeq2_pseudobulk_by_celltype",
    top_n = 50,
    min_donors_per_group = 2
  )
}, future.seed = TRUE)
#be careful here there is an eerror , te code mkes the plots though
names(all_results) <- names(pseudobulk_rna_results)
names(all_results) <-ct_vector
# optional: clean up
plan(sequential)

 
# keep only successful ones
all_results <- all_results[!sapply(all_results, is.null)]
length(all_results)




##############
#GSEA
############


# install if needed:
# install.packages(c("fgsea", "msigdbr", "dplyr", "tibble", "purrr", "ggplot2", "stringr"))
# BiocManager::install("fgsea")  # if needed

library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(stringr)

# Gene sets: Hallmark
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

pathways_h <- split(hallmark$gene_symbol, hallmark$gs_name)

# Gene sets: Reactome (MSigDB C2:CP:REACTOME)
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, gene_symbol)

pathways_react <- split(reactome$gene_symbol, reactome$gs_name)

# OptionA : use only Hallmark first (2 B faster)
pathways <- pathways_h
#pathways <- c(pathways_h, pathways_react)


# Helper: build a ranked vector + run fgsea
#Best practice: rank by Wald statistic from DESeq2 (directional + stable).deduplicate genes if needed and drop NAs.

make_rank <- function(res_df, rank_col = "stat") {
  # res_df: data.frame with rownames = genes
  stopifnot(rank_col %in% colnames(res_df))
  
  ranks <- res_df[[rank_col]]
  names(ranks) <- rownames(res_df)
  
  # remove NA / infinite
  ranks <- ranks[is.finite(ranks)]
  
  # if duplicated gene names exist, keep the one with largest |rank|
  if (any(duplicated(names(ranks)))) {
    ranks <- tibble(gene = names(ranks), rank = as.numeric(ranks)) %>%
      group_by(gene) %>%
      slice_max(order_by = abs(rank), n = 1, with_ties = FALSE) %>%
      ungroup()
    ranks <- ranks$rank
    names(ranks) <- ranks$gene
  }
  
  # sort decreasing for fgsea
  sort(ranks, decreasing = TRUE)
}
# Run GSEA across all cell types in all_results
run_fgsea_one <- function(ranks, pathways, minSize = 15, maxSize = 500) {
  fgseaMultilevel(
    pathways = pathways,
    stats    = ranks,
    minSize  = minSize,
    maxSize  = maxSize
  ) %>%
    as_tibble() %>%
    arrange(padj, desc(abs(NES)))
}

gsea_results <- imap(all_results, function(obj, ct) {
  res_df <- as.data.frame(obj$res)
  
  # DESeq2 results object -> data.frame keeps rownames as genes
  # Ensure 'stat' exists (it should)
  if (!("stat" %in% colnames(res_df))) {
    message("Skipping ", ct, ": no 'stat' column found.")
    return(NULL)
  }
  
  ranks <- make_rank(res_df, rank_col = "stat")
  
  fg <- run_fgsea_one(ranks, pathways = pathways, minSize = 15, maxSize = 500) %>%
    mutate(celltype = ct)
  
  fg
})

gsea_results <- gsea_results[!sapply(gsea_results, is.null)]
gsea_all <- bind_rows(gsea_results)

# Save long table
#dir.create("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/GSEA_results", showWarnings = FALSE, recursive = TRUE)
saveRDS(gsea_all, "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/GSEA_results/GSEA_all_celltypes_Hallmark_Reactome.csv")

# Make a pathway × cell type dot-plot (NES + FDR)

# Choose top pathways overall (by best padj across any cell type)
top_pathways <- gsea_all %>%
  filter(!is.na(padj)) %>%
  group_by(pathway) %>%
  summarise(best_padj = min(padj), .groups = "drop") %>%
  arrange(best_padj) %>%
  slice_head(n = 30) %>%
  pull(pathway)

plot_df <- gsea_all %>%
  filter(pathway %in% top_pathways) %>%
  mutate(
    neglog10_fdr = -log10(padj + 1e-300),
    pathway = factor(pathway, levels = rev(top_pathways))
  )
#Make a pathway × cell type dot-plot (NES + FDR)
nesfdr = ggplot(plot_df, aes(x = celltype, y = pathway)) +
  geom_point(aes(size = neglog10_fdr, color = NES)) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(linewidth = 0.2)
  ) +
  labs(
    title = "GSEA across cell types (Hallmark + Reactome)",
    x = "Cell type",
    y = NULL,
    size = "-log10(FDR)",
    color = "NES"
  )
ggsave(
  filename = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/GSEA_results/GSEA_dotplot_NES_FDR_hallmarkonly.png",
  plot = nesfdr,
  width = 12, height = 8, units = "in", dpi = 300
)
#nk cell investigation 

# Idents(seur_obj) <- "celltype_simplified"
# DefaultAssay(seur_obj) <- "RNA"
# DotPlot(
#   seur_obj, group.by = "sctype_classification_man",
#   features = c(
#     # NK
#     "NKG7","GNLY","PRF1","GZMB","KLRD1","FCGR3A","TRDC"
#   )
# ) + RotatedAxis()
# 
# DotPlot(
#   seur_obj, group.by = "sctype_classification_man",
#   features = c(
#     # T cell / TCR
#     "TRAC","TRBC1","CD3D","CD3E","LCK","CD8A","CD8B"
#   )
# ) + RotatedAxis()
# 
# 
# "sctype_classification_man"
# # NK subset
# nk <- subset(seur_obj, subset = sctype_classification_man == "Natural killer  cells")
# #
# DefaultAssay(nk) <- "RNA"
# nk_counts <- GetAssayData(nk, layer = "counts")
# 
# tcr_genes <- c("TRAC","TRBC1","CD3D","CD3E")
# tcr_genes <- tcr_genes[tcr_genes %in% rownames(nk_counts)]
# 
# tcr_frac <- sapply(tcr_genes, function(g) mean(nk_counts[g, ] > 0))
# tcr_frac
# 
# 
# FeatureScatter(nk, feature1 = "NKG7", feature2 = "TRAC")
# FeatureScatter(nk, feature1 = "GNLY", feature2 = "TRAC")
# DotPlot(seur_obj, features = tcell_rna_markers, group.by = "sctype_classification_man",
#         cols = c("lightgrey", "blue"), assay = "RNA") +
#   RotatedAxis() +
#   ggtitle("Canonical T cell marker expression")
# 
# VlnPlot(seur_obj, features=c("CD3","CD56","CD16"), group.by="sctype_classification_man", pt.size=0)
# 
#  
# 
# nk[["RNA"]] <- JoinLayers(
#   object = nk[["RNA"]],
#   layers = "counts"
# )
# 
# 
# 
# nk <- subset(seur_obj, subset = celltype_simplified == "NK cells")
# DefaultAssay(nk) <- "RNA"
# nk_counts <- GetAssayData(nk, layer = "counts")
# 
# t_genes <- c("TRAC","TRBC1","CD3D","CD3E","LCK")
# t_genes <- t_genes[t_genes %in% rownames(nk_counts)]
# 
# t_umis <- Matrix::colSums(nk_counts[t_genes, , drop=FALSE])
# 
# summary(t_umis)
# mean(t_umis > 0)
# mean(t_umis >= 2)
# mean(t_umis >= 5)
# 
# # per-gene quick check (how many UMIs when present)
# for (g in t_genes) {
#   v <- nk_counts[g, ]
#   cat(g, ": median(nonzero)=", median(v[v>0]), "  max=", max(v), "\n")
# }
# 
# # strict: no TCR/CD3 UMIs at all
# nk_pure_cells <- colnames(nk)[t_umis == 0]
# nk_pure <- subset(nk, cells = nk_pure_cells)
# 
# dim(nk_pure)
# 
# 
# nk_pure_counts <- GetAssayData(nk_pure, layer="counts")
# 
# tcr_genes <- c("TRAC","TRBC1","CD3D","CD3E","LCK")
# tcr_genes <- tcr_genes[tcr_genes %in% rownames(nk_pure_counts)]
# 
# tcr_frac_pure <- sapply(tcr_genes, function(g) mean(nk_pure_counts[g, ] > 0))
# tcr_frac_pure
# 
# 
# table(nk_pure$validated_TB_status)
# length(unique(nk_pure$NCR.ID))
# 
# cells_per_donor <- table(nk_pure$NCR.ID)
# summary(as.numeric(cells_per_donor))
# 
# # How many donors have at least N NK cells?
# sum(cells_per_donor >= 20)
# sum(cells_per_donor >= 50)
# sum(cells_per_donor >= 100)
# 
# #Make a “NK_pure” label in the full Seurat object (so your existing pseudobulk function still works)
# # 1) start everyone as their current simplified type
# seur_obj$celltype_simplified2 <- seur_obj$celltype_simplified
# 
# # 2) mark the pure NK cells
# seur_obj$celltype_simplified2[colnames(seur_obj) %in% colnames(nk_pure)] <- "NK_pure"
# 
# # optional: keep the remaining NK cells as "NK_mixed"
# nk_all_cells <- WhichCells(seur_obj, expression = celltype_simplified == "NK cells")
# nk_pure_cells <- colnames(nk_pure)
# nk_mixed_cells <- setdiff(nk_all_cells, nk_pure_cells)
# seur_obj$celltype_simplified2[nk_mixed_cells] <- "NK_mixed"
# 
# DotPlot(
#   seur_obj,
#   features = c("NKG7","GNLY","PRF1","GZMB","KLRD1","FCGR3A",
#                "TRAC","TRBC1","CD3D","CD3E","LCK",
#                "CD8A","CD8B"), group.by = "sctype_classification_man"
# ) + RotatedAxis()
