suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})
set.seed(626)
 
 
library(dplyr)
 
library(ggalluvial)
library(patchwork)
rds_path <-"/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/ScType_multiDB_out_res12_pca15_oscar/raw_merge_all_batches_harm_annotated_all_res12_pca20_noann_Seurat_ScType_6DB_oscar.rds"

 #rds_path <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/ScType_multiDB_out_noTCRres12_pca15/raw_merge_all_batches_harm_annotated_all_res12_pca15_noann_Seurat_ScType_7DB.rds"

obj <- readRDS(rds_path)
 
meta <- obj@meta.data
  
  
  # 

 

DefaultAssay(obj) <- "RNA"
Reductions(obj)
colnames(obj@meta.data)[1:45]
nrow(obj)

 
donor_col    <- "NCR.ID"
status_col   <- "validated_TB_status"
batch_col    <- "group"   
table(obj@meta.data[[celltype_col]]) %>% sort(decreasing = TRUE) %>% head(30)

#subset to nk cells
#nk <- subset(obj, subset = sctype_cd8_nk_tcr_neg == "Natural killer  cells")
nk <- subset(obj, subset = sctype_default == "Natural killer  cells")
#nk <- subset(obj, subset = sctype_no_isg == "Natural killer  cells")

 
print(ElbowPlot(nk, ndims = 60))
# Use harmony embedding if present
reduction_use <- if ("harmony" %in% Reductions(nk)) "harmony" else "pca"
reduction_use

# Pick dims (you used pca35; 1:35 is a reasonable default)
dims_use <- 1:35

nk <- FindNeighbors(nk, reduction = reduction_use, dims = dims_use)
nk <- FindClusters(nk, resolution = 0.6)   # try 0.3–1.0 if you want
nk <- RunUMAP(nk, reduction = reduction_use, dims = dims_use, reduction.name = "umap_nk")


nk <- subset(obj, subset = sctype_default == "Natural killer  cells")
DefaultAssay(nk) <- "RNA"

nk <- FindVariableFeatures(nk, nfeatures = 2000)
nk <- ScaleData(nk)
nk <- RunPCA(nk, npcs = 60)

nk <- harmony::RunHarmony(nk, group.by.vars = "batch", reduction = "pca")
nk<- RunHarmony(
  object = nk,
  group.by.vars = "batch",
  dims.use = 1:20,
  assay.use = "RNA",
  plot_convergence = TRUE
)


nk <- FindNeighbors(nk, reduction = "harmony", dims = 1:20)
nk <- FindClusters(nk, resolution = 0.6)
nk <- RunUMAP(nk, reduction = "harmony", dims = 1:20, reduction.name = "umap_nk")


DimPlot(nk, reduction = "umap_nk", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("NK-only reclustering (Harmony space)")

DimPlot(nk, reduction = "umap_nk", group.by = batch_col) +
  ggtitle("NK UMAP colored by batch")
 

nk1 <- c(
  "CD160",
  "CTSD",
  "CCL4",
  "ADGRG1",
  "CD38",
  "CD247",
  "CHST2",
  "CX3CR1",
  "KLRB1",
  "LAIR2",
  "IGFBP7",
  "AKR1C3",
  "FGFBP2",
  "MYOM2",
  "CLIC3",
  "GZMB",
  "PRF1",
  "FCER1G",
  "NKG7",
  "SPON2"
)
nk2 <- c(
  "LTB",
  "FOS",
  "IL2RB",
  "IFITM3",
  "COTL1",
  "IL7R",
  "PIK3R1",
  "AREG",
  "ZFP36L2",
  "DUSP2",
  "CD44",
  "SELL",
  "GPR183",
  "CMC1",
  "KLRC1",
  "TCF7",
  "TPT1",
  "XCL2",
  "XCL1",
  "GZMK"
)
nk3 <- c(
  "HBA1",
  "CD3G",
  "CD3D",
  "DSTN",
  "ZBTB38",
  "CD2",
  "PPDPF",
  "PRDM1",
  "LGALS1",
  "S100A4",
  "ITGB1",
  "IL32",
  "VIM",
  "LINC01871",
  "S100A6",
  "PTMS",
  "GZMH",
  "CD3E",
  "CCL5",
  "KLRC2"
)

sig_A <- intersect(nk1, rownames(nk))
sig_B <- intersect(nk2, rownames(nk))
sig_C <- intersect(nk3, rownames(nk))



length(sig_A); length(sig_B); length(sig_C)


nk <- AddModuleScore(nk, features = list(sig_A), name = "NKsigA_")
nk <- AddModuleScore(nk, features = list(sig_B), name = "NKsigB_")
nk <- AddModuleScore(nk, features = list(sig_C), name = "NKsigC_")

head(nk@meta.data[, c("NKsigA_1","NKsigB_1","NKsigC_1")])


FeaturePlot(nk, reduction = "umap_nk", features = c("NKsigA_1","NKsigB_1","NKsigC_1"), ncol = 3)
features_use <- list(
  NKsigA = sig_A,
  NKsigB = sig_B,
  NKsigC = sig_C
)

DotPlot(nk, features = features_use) + RotatedAxis()

DotPlot(nk, features = c("NKsigA_1","NKsigB_1","NKsigC_1"),
        group.by = "sctype_cd8_nk") + RotatedAxis()
scores <- nk@meta.data[, c("NKsigA_1","NKsigB_1","NKsigC_1")]
colnames(scores) <- c("A","B","C")

top <- apply(scores, 1, function(x) names(which.max(x)))
sorted2 <- t(apply(scores, 1, function(x) sort(x, decreasing = TRUE)[1:2]))
delta <- sorted2[,1] - sorted2[,2]

margin <- 0.15   # tweak: 0.1–0.25 typical for AddModuleScore
nk$nk_subset3 <- ifelse(delta >= margin, top, "Ambiguous")
nk$nk_subset3 <- factor(nk$nk_subset3, levels = c("A","B","C","Ambiguous"))

table(nk$nk_subset3)
DimPlot(nk, reduction = "umap_nk", group.by = "nk_subset3") + ggtitle("NK subset assignment (Default)")




meta <- nk@meta.data %>%
  mutate(
    donor  = .data[[donor_col]],
    status = .data[[status_col]]
  ) %>%
  filter(!is.na(donor), !is.na(status))

 
donor_counts <- meta %>%
  filter(nk_subset3 != "Ambiguous") %>%
  count(donor, status, Gender,group, nk_subset3, name = "n_cells") %>%  # <-- add sex
  group_by(donor) %>%
  mutate(nk_total = sum(n_cells),
         prop = n_cells / nk_total) %>%
  ungroup()

 
library(ggplot2)
library(ggpubr)

ggplot(donor_counts, aes(x = status, y = prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Gender), width = 0.15, height = 0, alpha = 0.7) +
  facet_wrap(~ nk_subset3) +
  theme_bw(base_size = 12) +
  labs(title = "Donor-level NK subset proportions (CD8/NK(TCR-))", x = NULL, y = "Proportion") +
  stat_compare_means(method = "wilcox.test", label = "p.format")

# Horizontal bars of group means (your earlier style)
group_means <- donor_counts %>%
  group_by(status, nk_subset3) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

ggplot(group_means, aes(x = mean_prop, y = nk_subset3, fill = status)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  theme_bw(base_size = 12) +
  labs(title = "Mean NK subset proportion by group (donor-level)", x = "Mean proportion", y = NULL)
adt_markers <- c(
  "Hu.CD27",
  "Hu.CD2",
  "Hu.CD56",
  "Hu.CD314",
  "Hu.CD117",
  "Hu.CD96",
  "Hu.CD54",
  "Hu.CD45RB",
  "Hu.CD335",
  "Hu.CD44",
  "Hu.CD38_HIT2","Hu.CD122","Hu.CD16","Hu.CD161","Hu.CD337"
)

DefaultAssay(nk) <- "ADT"   # or whatever your protein assay is called

DotPlot(
  nk,
  features = adt_markers,
  group.by = "nk_subset3"   # change if you want another x-axis grouping
) + RotatedAxis()



DotPlot(
  nk,
  features = "ANXA2R",
  group.by = "nk_subset3"   # change if you want another x-axis grouping
) + RotatedAxis()




DotPlot(
  obj,
  features = "ANXA2R",
  group.by = "sctype_default","sctype_cd8_nk_tcr_neg", "sctype_cd8_nk"
) +
  RotatedAxis()

FeaturePlot(
  obj,
  features = "ANXA2R",
  reduction = "rna.umap",
  cols = c("lightgrey", "red")
)


p1 <- DotPlot(
  obj,
  features = "ANXA2R",
  group.by = "sctype_default"
) + RotatedAxis() + ggtitle("Default annotation")

p2 <- DotPlot(
  obj,
  features = "ANXA2R",
  group.by = "sctype_cd8_nk_tcr_neg"
) + RotatedAxis() + ggtitle("CD8/NK (TCR−)")

p3 <- DotPlot(
  obj,
  features = "ANXA2R",
  group.by = "sctype_cd8_nk"
) + RotatedAxis() + ggtitle("CD8/NK")

 
p1 | p2 | p3
p1 | p2 


p_case <- DotPlot(
  subset(obj, subset = validated_TB_status == "2weekCase"),
  features = "ANXA2R",
  group.by = "sctype_default"
) +
  RotatedAxis() +
  ggtitle("Default annotation— Cases")

p_ctrl <- DotPlot(
  subset(obj, subset = validated_TB_status == "ctrl"),
  features = "ANXA2R",
  group.by = "sctype_default"
) +
  RotatedAxis() +
  ggtitle("Default annotation — Controls")


p_case | p_ctrl


 

p_case_vln <- VlnPlot(
  subset(obj, subset = validated_TB_status == "2weekCase"),
  features = "ANXA2R",
  group.by = "sctype_default",
  pt.size = 0
) + ggtitle("ANXA2R — Cases") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ctrl_vln <- VlnPlot(
  subset(obj, subset = validated_TB_status == "ctrl"),
  features = "ANXA2R",
  group.by = "sctype_default",
  pt.size = 0
) + ggtitle("ANXA2R — Controls") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_case_vln | p_ctrl_vln

Idents(obj) <- "sctype_default"

p_split <- VlnPlot(
  obj_filt,
  features = "ANXA2R",
  group.by = "sctype_default",
  split.by = "validated_TB_status",
  pt.size = 0
) + ggtitle("ANXA2R — split by TB status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 p_split
ggsave(
  filename = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/p_split_violin.png",
  plot = p_split,
  width = 13,
  height = 5,
  dpi = 300
)

p_split <- VlnPlot(
  obj,
  features = "ANXA2R",
  group.by = "sctype_default",
  split.by = "validated_TB_status",
  pt.size = 0,
  fill.by = "split"
) +
  scale_fill_manual(values = c("2weekCase" = "#E76F51", "ctrl" = "#2A9D8F")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_split


 