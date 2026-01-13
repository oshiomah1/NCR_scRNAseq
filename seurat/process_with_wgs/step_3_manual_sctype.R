lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

library(Seurat)
library(ggplot2)
library(scales)
library(dplyr)
library(openxlsx)
library(HGNChelper)
set.seed(411)
library(patchwork)


#seur_obj <- readRDS("quobyte/from-lssc0/projects/NCR_scRNAseq/results/seurat/harmony_22jan_.rds") #reg processing + wnn
seur_obj <- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/raw_merge_all_batches_harm_annotated_all_res12_pca35_noann.rds") #reg processing + wnn

print("read file")
head(seur_obj)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# DB file
#db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx";
#db_ <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/ScTypeDB_immune_only.xlsx"
#db_ <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/no_isg.xlsx"
#db_ <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/no_isg_no_mega.xlsx" #man5
db_ <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/no_isg_nktlikenegedit.xlsx" # man 3
#db_ <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/no_isg_nktlikepozedit.xlsx" # man 4
 

db_ <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/no_isg_nktlikenegedit_nomega.xlsx" #man 6 is 
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)



# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seur_obj[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(seur_obj[["RNA"]]$scale.data) else as.matrix(seur_obj[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. seur_obj[["RNA"]]$counts or seur_obj[["RNA"]]@counts, with scaled set to FALSE.

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(seur_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seur_obj@meta.data[seur_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seur_obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])
 
seur_obj@meta.data$sctype_classification_man = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  seur_obj@meta.data$sctype_classification_man[seur_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


umap = DimPlot(seur_obj, reduction = "rna.umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification_man2')        
DimPlot(seur_obj, reduction = "rna.umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification_man')        
library(ggplot2)
DimPlot(seur_obj, reduction = "rna.umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification_man6')       
 

ggsave(
  filename = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/rna_umap_sctype_man_res1.png",
  plot = umap,
  width = 14, height = 7, units = "in", dpi = 300
)

saveRDS(seur_obj, file = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/raw_merge_all_batches_harm_sc_annotated_res12_pca35_multianntest.rds")

print("file_saved")

# prepare edges
cL_resutls <- cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 <- sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss <- c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes <- rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db <- openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]


top4 <- cL_resutls %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(scores, n = 4, with_ties = FALSE) %>%
  dplyr::ungroup()
library(ggplot2)


pdf("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/ScType_top4_per_cluster_pcs15_res1_notcr.pdf", width = 7, height = 5)
 

for (cl in sort(unique(top4$cluster))) {
  
  df_cl <- dplyr::filter(top4, cluster == cl)
  
  p <- ggplot(df_cl, aes(x = reorder(type, scores), y = scores, fill = type)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Cluster", cl, "- Top 4 ScType assignments"),
      x = NULL,
      y = "ScType score"
    )
  
  print(p)
}

dev.off()

markers <-read.csv("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/Markers.csv", na.strings = c("", "NA"))

tcell_rna_markers <- na.omit(markers[markers$Cell.type == "T-cells", "RNA_name"])
bcell_rna_markers <- na.omit(markers[markers$Cell.type == "B-cell", "RNA_name"])
monocytes_rna_markers <- na.omit(markers[markers$Cell.type == "Monocytes", "RNA_name"])
nkcell_rna_markers <- na.omit(markers[markers$Cell.type == "NK-cells", "RNA_name"])
dendritic_cell_rna_markers <- na.omit(markers[markers$Cell.type == "Dendritic-Cell", "RNA_name"])
adt_markers <- na.omit(markers$ADT_name)

DotPlot(
  seur_obj,
  features = c("CD8", "CD56", "CD2", "CD16", "CD94", "CD3D", "CD3E", "CD3G", "CD3Z",
               "NKp46", "CD11b", "CD161", "CD314", "CD69", "NKG7", "CD122", "NKG2D",
               "GZMB", "GZMA", "GZMM", "GNLY", "COX6A2", "ZMAT4", "KIR2DL4"), group.by = "sctype_classification_man2"
) + RotatedAxis()
PECAM1,CD34,KDR,CDH5,PROM1,PDPN,TEK,FLT1,VCAM1,PTPRC,VWF,ENG,MCAM,ICAM1,FLT4,SELE
DotPlot(
  seur_obj,
  features = c("PECAM1", "CD34", "KDR", "CDH5", "PROM1", "PDPN", "TEK", "FLT1",
               "VCAM1", "PTPRC", "VWF", "ENG", "MCAM", "ICAM1", "FLT4", "SELE"), group.by = "sctype_classification_man2"
) + RotatedAxis()
DotPlot(seur_obj, features = tcell_rna_markers, group.by = "sctype_classification_man6",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical T cell marker expression")
DotPlot(seur_obj, features = nkcell_rna_markers, group.by = "sctype_classification_man6",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical NK cell marker expression")
DotPlot(seur_obj, features = bcell_rna_markers, group.by = "sctype_classification_man6",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical B cell marker expression")

DotPlot(seur_obj, features = monocytes_rna_markers, group.by = "sctype_classification_man2",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical Monocyte marker expression")
DotPlot(seur_obj, features = adt_markers, group.by = "sctype_classification_man2",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("adt marker expression")
DotPlot(seur_obj, features = adt_markers, group.by = "sctype_classification_man2",
        cols = c("lightgrey", "blue"), assay = "ADT") +
  RotatedAxis() +
  ggtitle("Canonical Cell surface protein marker expression")

# Aggregate cell counts (no age/gender yet)
celltype_props <- seur_obj@meta.data %>%
  dplyr::filter(!is.na(NCR.ID)) %>%
  dplyr::count(NCR.ID, validated_TB_status, sctype_classification_man, name = "n_cells") %>%
  dplyr::group_by(NCR.ID) %>%
  dplyr::mutate(
    total_cells = sum(n_cells),
    prop = n_cells / total_cells
  ) %>%
  dplyr::ungroup()

# Next we extract donor-level metadata (age, gender, etc.)

donor_meta <- seur_obj@meta.data %>%
  dplyr::filter(!is.na(NCR.ID)) %>%
  dplyr::select(
    NCR.ID,
    Age.,
    Gender
  ) %>%
  dplyr::distinct()

#Join donor metadata onto proportions
celltype_props <- celltype_props %>%
  dplyr::left_join(donor_meta, by = "NCR.ID")

#rea in abx data
abx<- read.csv("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/extra_metadata.tsv") %>%
  rename(
    days_antibiotics = lab_date_2_trtment_strt
  ) %>%
  mutate(
    validated_case_status = factor(validated_case_status),
    days_antibiotics = as.numeric(days_antibiotics)
  )

celltype_props <- celltype_props %>%
  left_join(abx, by = "NCR.ID")

##horizontal bar plot

summary_df <- celltype_props %>%
  group_by(validated_TB_status, sctype_classification_man) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percent = mean_prop * 100
  )

celltype_order <- summary_df %>%
  group_by(sctype_classification_man) %>%
  summarise(overall_mean = mean(percent)) %>%
  arrange(overall_mean) %>%
  pull(sctype_classification_man)

summary_df$sctype_classification_man <-
  factor(summary_df$sctype_classification_man, levels = celltype_order)

#cell type prportions plus inset!

p_main <- ggplot(
  summary_df,
  aes(
    x = percent,
    y = sctype_classification_man,
    fill = validated_TB_status
  )
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(
    values = c(
      "2weekCase" = "#FFC20A",
      "ctrl"      = "#0C7BDC"
    ),
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
  labs(
    x = "Cell type proportion (%)",
    y = NULL
  )

#select the cell types with low proportions
bottom_types <- head(levels(summary_df$sctype_classification_man), 3)

summary_df_small <- summary_df %>%
  dplyr::filter(sctype_classification_man %in% bottom_types)

#plot  the cell types with low proportions as an inset
p_inset <- ggplot(
  summary_df_small,
  aes(
    x = percent,
    y = sctype_classification_man,
    fill = validated_TB_status
  )
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(
    values = c(
      "2weekCase" = "#FFC20A",
      "ctrl"      = "#0C7BDC"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 0.8),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8),
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
  labs(
    x = "Zoomed (%)",
    y = NULL
  )

library(patchwork)

#use patchwork to stutch inset and main
barplots_celltypes_inset <- p_main + inset_element(
  p_inset,
  left   = 0.45,
  bottom = 0.08,
  right  = 0.95,
  top    = 0.50
)


ggsave(
  filename = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/cell_type_proportions/barplots_celltypes_inset.png",
  plot     = barplots_celltypes_inset,
  width    = 18,
  height   = 6,
  units    = "in",
  dpi      = 300
)

md <- seur_obj@meta.data

# sanity: what VDJ-related columns exist?
grep("cdr3|clonotype|chain|vdj|tcr|bcr", colnames(md), ignore.case = TRUE, value = TRUE)

# define "has TCR" and "has BCR" using whichever fields you actually have
has_tcr <- rep(FALSE, nrow(md))
has_bcr <- rep(FALSE, nrow(md))

if ("t_clonotype_id" %in% colnames(md)) has_tcr <- has_tcr | !is.na(md$t_clonotype_id) & md$t_clonotype_id != ""
if ("t_cdr3s_aa"     %in% colnames(md)) has_tcr <- has_tcr | !is.na(md$t_cdr3s_aa)     & md$t_cdr3s_aa != ""

if ("b_clonotype_id" %in% colnames(md)) has_bcr <- has_bcr | !is.na(md$b_clonotype_id) & md$b_clonotype_id != ""
if ("b_cdr3s_aa"     %in% colnames(md)) has_bcr <- has_bcr | !is.na(md$b_cdr3s_aa)     & md$b_cdr3s_aa != ""

md$has_tcr <- has_tcr
md$has_bcr <- has_bcr

# overall counts
table(md$has_tcr, md$has_bcr)

# suspicious multiplets: both TCR and BCR
mean(md$has_tcr & md$has_bcr)

# now: enrichment by  ScType labels
prop_by_type <- with(md, tapply(has_tcr, sctype_classification_man6, mean))
sort(prop_by_type, decreasing = TRUE)

prop_by_type_bcr <- with(md, tapply(has_bcr, sctype_classification_man6, mean))
sort(prop_by_type_bcr, decreasing = TRUE)





# contingency table
print(table(md$sctype_classification_man6, md$has_tcr))

# fraction TCR+ per label (tidy)
library(dplyr)
md %>%
  group_by(sctype_classification_man6) %>%
  summarise(
    n = n(),
    frac_tcr = mean(has_tcr),
    frac_bcr = mean(has_bcr),
    frac_both = mean(has_tcr & has_bcr)
  ) %>%
  arrange(desc(frac_tcr))

tcr_vec <- md$t_cdr3s_aa
tcr_vec <- tcr_vec[!is.na(tcr_vec) & tcr_vec != ""]

tcr_freq <- sort(table(tcr_vec), decreasing = TRUE)
top_10_t <- head(tcr_freq[tcr_freq >= 10], 10)


# # Build matrix: rows = cell types, cols = clusters
# mat <- cL_resutls %>%
#   dplyr::select(cluster, type, scores) %>%
#   tidyr::pivot_wider(names_from = cluster, values_from = scores, values_fill = 0) %>%
#   tibble::column_to_rownames("type") %>%
#   as.matrix()
# 
# pheatmap(
#   mat,
#   scale = "row",
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   fontsize_row = 7,
#   fontsize_col = 9
# )
