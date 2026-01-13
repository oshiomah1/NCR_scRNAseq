suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})
set.seed(626)
ee
#install.packages("ggalluvial")  # if needed
library(dplyr)
library(ggplot2)
library(ggalluvial)

rds_path <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/raw_merge_all_batches_harm_annotated_all_res12_pca35_noann_Seurat_ScType_7DB.rds"
#rds_path <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/ScType_multiDB_out_noTCRres12_pca15/raw_merge_all_batches_harm_annotated_all_res12_pca15_noann_Seurat_ScType_7DB.rds"
obj <- readRDS(rds_path)

meta <- obj@meta.data
#################################
#########VDJ CHECK###############
#################################



# sanity: what VDJ-related columns exist?
grep("cdr3|clonotype|chain|vdj|tcr|bcr", colnames(meta), ignore.case = TRUE, value = TRUE)

# define "has TCR" and "has BCR" using whichever fields you actually have
has_tcr <- rep(FALSE, nrow(meta))
has_bcr <- rep(FALSE, nrow(meta))

if ("t_clonotype_id" %in% colnames(meta)) has_tcr <- has_tcr | !is.na(meta$t_clonotype_id) & meta$t_clonotype_id != ""
if ("t_cdr3s_aa"     %in% colnames(meta)) has_tcr <- has_tcr | !is.na(meta$t_cdr3s_aa)     & meta$t_cdr3s_aa != ""

if ("b_clonotype_id" %in% colnames(meta)) has_bcr <- has_bcr | !is.na(meta$b_clonotype_id) & meta$b_clonotype_id != ""
if ("b_cdr3s_aa"     %in% colnames(meta)) has_bcr <- has_bcr | !is.na(meta$b_cdr3s_aa)     & meta$b_cdr3s_aa != ""

meta$has_tcr <- has_tcr
meta$has_bcr <- has_bcr

# overall counts
table(meta$has_tcr, meta$has_bcr)

# suspicious multiplets: both TCR and BCR
mean(meta$has_tcr & meta$has_bcr)

# now: enrichment by  ScType labels
prop_by_type <- with(meta, tapply(has_tcr, sctype_no_isg , mean))
sort(prop_by_type, decreasing = TRUE)

prop_by_type_bcr <- with(meta, tapply(has_bcr, sctype_baseline, mean))
sort(prop_by_type_bcr, decreasing = TRUE)


#################################
#########VDJ CHECK###############
#################################


cols <- c(
  "sctype_baseline",
  "sctype_nktlikenegedit",
  "sctype_no_isg"
)

lump_n <- 12  # keep top 12 labels per column

meta2 <- meta %>%
  select(all_of(cols)) %>%
  mutate(across(everything(), \(x) {
    x <- as.character(x)
    keep <- names(sort(table(x), decreasing = TRUE))[1:min(lump_n, length(unique(x)))]
    ifelse(x %in% keep, x, "Other")
  }))

ggplot(meta2,
       aes(axis1 = sctype_baseline,
           axis2 = sctype_nktlikenegedit,
           axis3 = sctype_no_isg
          )) +
  geom_alluvium(aes(fill = sctype_baseline), alpha = 0.8) +
  geom_stratum(color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = cols, expand = c(.02, .02)) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Baseline")


library(ggplot2)
library(ggalluvial)

p <- ggplot(meta2,
            aes(axis1 = sctype_baseline,
                axis2 = sctype_nktlikenegedit,
                axis3 = sctype_no_isg)) +
  geom_alluvium(aes(fill = sctype_baseline), alpha = 0.8) +
  geom_stratum(color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("sctype_baseline","sctype_nktlikenegedit","sctype_no_isg"),
                   expand = c(.02, .02)) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Baseline")

# Save it much wider
ggsave("alluvial_wide.png", plot = p, width = 16, height = 6, units = "in", dpi = 300)
# or PDF (nice for text)
ggsave("alluvial_wide.pdf", plot = p, width = 16, height = 6, units = "in")

  
  #################################
  #########VDJ CHECK###############
  #################################
  







DefaultAssay(obj) <- "RNA"
Reductions(obj)
colnames(obj@meta.data)[1:45]
nrow(obj)

 
donor_col    <- "NCR.ID"
status_col   <- "validated_TB_status"
batch_col    <- "group"   
table(obj@meta.data[[celltype_col]]) %>% sort(decreasing = TRUE) %>% head(30)

#subset to nk cells
nk <- subset(obj, subset = sctype_nktlikenegedit == "Natural killer  cells")
#nk <- subset(obj, subset = sctype_baseline == "Natural killer  cells")
#nk <- subset(obj, subset = sctype_no_isg == "Natural killer  cells")

 

# Use harmony embedding if present
reduction_use <- if ("harmony" %in% Reductions(nk)) "harmony" else "pca"
reduction_use

# Pick dims (you used pca35; 1:35 is a reasonable default)
dims_use <- 1:35

nk <- FindNeighbors(nk, reduction = reduction_use, dims = dims_use)
nk <- FindClusters(nk, resolution = 0.6)   # try 0.3–1.0 if you want
nk <- RunUMAP(nk, reduction = reduction_use, dims = dims_use, reduction.name = "umap_nk")
DimPlot(nk, reduction = "umap_nk", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("NK-only reclustering (Harmony space)")

DimPlot(nk, reduction = "umap_nk", group.by = batch_col) +
  ggtitle("NK UMAP colored by batch")
#nk2  <- c("IL2RB", "TPT1", "GPR183", "CD44", "TCF7",
              #   "IL7R", "KLRC1", "XCL1", "SELL", "CMC1",
               #  "XCL2", "GZMK", "CD2")
#nk3  <- c("CD52", "S100A4", "LGALS1", "PTMS", "LINC01871",
             #   "VIM", "S100A6", "ITGB1", "IL32", "CCL5",
             #   "GZMH", "CD3E", "KLRC2")
#nk1  <- c("PLAC8", "CD38", "CCL4", "CHST2", "FCER1G", "IGFBP7", "AKR1C3", "SPON2")


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
        group.by = "sctype_nktlikenegedit") + RotatedAxis()
scores <- nk@meta.data[, c("NKsigA_1","NKsigB_1","NKsigC_1")]
colnames(scores) <- c("A","B","C")

top <- apply(scores, 1, function(x) names(which.max(x)))
sorted2 <- t(apply(scores, 1, function(x) sort(x, decreasing = TRUE)[1:2]))
delta <- sorted2[,1] - sorted2[,2]

margin <- 0.15   # tweak: 0.1–0.25 typical for AddModuleScore
nk$nk_subset3 <- ifelse(delta >= margin, top, "Ambiguous")
nk$nk_subset3 <- factor(nk$nk_subset3, levels = c("A","B","C","Ambiguous"))

table(nk$nk_subset3)
DimPlot(nk, reduction = "umap_nk", group.by = "nk_subset3") + ggtitle("NK subset assignment (A/B/C)")





meta <- obj@meta.data %>%
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
  labs(title = "Donor-level NK subset proportions (within NK)", x = NULL, y = "Proportion") +
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






sctype_baseline




 