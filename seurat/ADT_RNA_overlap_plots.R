library(Seurat)
#quobyte/from-lssc0/
library(ggplot2)
library(scales)
library(dplyr)
library(openxlsx)
library(HGNChelper)
set.seed(411)


library(patchwork)


#seur_obj <- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/raw_merge_all_batches_annotated.rds") #reg processing + wnn
seur_obj <- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/raw_merge_all_batches_harm_annotated_all.rds") #reg processing + wnn
seur_obj = obj
markers <-read.csv("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/Markers.csv", na.strings = c("", "NA"))
 
tcell_rna_markers <- na.omit(markers[markers$Cell.type == "T-cells", "RNA_name"])
bcell_rna_markers <- na.omit(markers[markers$Cell.type == "B-cell", "RNA_name"])
monocytes_rna_markers <- na.omit(markers[markers$Cell.type == "Monocytes", "RNA_name"])
nkcell_rna_markers <- na.omit(markers[markers$Cell.type == "NK-cells", "RNA_name"])
dendritic_cell_rna_markers <- na.omit(markers[markers$Cell.type == "Dendritic-Cell", "RNA_name"])
adt_markers <- na.omit(markers$ADT_name)

tcell_rna_markers 
bcell_rna_markers
monocytes_rna_markers
nkcell_rna_markers
dendritic_cell_rna_markers

my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#5254a3", "#9c9ede", "#d6616b", "#b00b69", "#de9ed6"
)
adt_markers <- na.omit(markers$ADT_name)

adt_markers

DotPlot(seur_obj, features = tcell_rna_markers, group.by = "sctype_classification",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical T cell marker expression")
DotPlot(seur_obj, features = bcell_rna_markers, group.by = "sctype_classification",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical B cell marker expression")

DotPlot(seur_obj, features = monocytes_rna_markers, group.by = "sctype_classification",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical Monocyte marker expression")

DotPlot(seur_obj, features = nkcell_rna_markers, group.by = "sctype_classification3",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical NK cell marker expression")

DotPlot(seur_obj, features = dendritic_cell_rna_markers, group.by = "sctype_classification",
        cols = c("lightgrey", "blue"), assay = "RNA") +
  RotatedAxis() +
  ggtitle("Canonical Dendritic cell marker expression")

DotPlot(seur_obj, features = adt_markers, group.by = "sctype_classification",
        cols = c("lightgrey", "blue"), assay = "ADT") +
  RotatedAxis() +
  ggtitle("Canonical Cell surface protein marker expression")

DotPlot(seur_obj, features = tcell_rna_markers, group.by = "sctype_classification", cols = c("lightgrey", "blue"), assay = "RNA") + RotatedAxis()  
DotPlot(seur_obj, features = bcell_rna_markers, group.by = "sctype_classification", cols = c("lightgrey", "blue"), assay = "RNA") + RotatedAxis()  
DotPlot(seur_obj, features = monocytes_rna_markers, group.by = "sctype_classification", cols = c("lightgrey", "blue"), assay = "RNA") + RotatedAxis()  
DotPlot(seur_obj, features = nkcell_rna_markers, group.by = "sctype_classification", cols = c("lightgrey", "blue"), assay = "RNA") + RotatedAxis()  
DotPlot(seur_obj, features = dendritic_cell_rna_markers, group.by = "sctype_classification", cols = c("lightgrey", "blue"), assay = "RNA") + RotatedAxis()  

DotPlot(seur_obj, features = adt_markers, group.by = "sctype_classification", assay ="ADT",cols = c("lightgrey", "blue")) + RotatedAxis()


DimPlot(
  object = seur_obj,
  reduction = use_red,
  group.by = "sctype_default",
  label = TRUE,             # place labels at centroids
  repel = TRUE,             # spread labels if many groups
  label.size = 3.5
) + ggtitle(paste("UMAP colored by sctype_classification (", use_red, ")", sep = "")) +
  theme_bw() + theme(panel.grid = element_blank())

DimPlot(seur_obj, reduction = 'rna.umap', group.by = "sctype_classification",  label = TRUE, repel = TRUE, label.size = 5.5) 


# 1. Pick RNA-only UMAP as the base
DimPlot(
  seur_obj,
  reduction = "rna.umap",
  group.by = "sctype_default",
  label = TRUE,
  repel = TRUE,
  label.size = 5.5
) + ggtitle("RNA UMAP colored by annotated cell types")

# 2. Overlay protein expression directly on the RNA UMAP
DefaultAssay(seur_obj) <- "ADT"

# Example:  
FeaturePlot(
  seur_obj,
  features = c("Hu.CD20-2H7", "Hu.CD19"),   # swap with your adt_markers
  reduction = "rna.umap",
  min.cutoff = "q05", max.cutoff = "q95",
  cols = c("lightgrey", "blue")
)

for (f in adt_markers) {
  print(
    FeaturePlot(
      seur_obj,
      features = f,
      reduction = "rna.umap",
      min.cutoff = "q05", max.cutoff = "q95",
      cols = c("lightgrey", "blue")
    ) + ggtitle(paste("ADT:", f, "on RNA UMAP"))
  )
}


pdf("ADT_markers_rnaUMAP_v2.pdf", width = 6, height = 5)

for (f in adt_markers) {
  p <- FeaturePlot(
    seur_obj,
    features = f,
    reduction = "rna.umap",
    min.cutoff = "q05", max.cutoff = "q95",
    cols = c("lightgrey", "blue")
  ) + ggtitle(paste("ADT:", f, "on RNA UMAP"))
  
  print(p)   # needed inside pdf loop
}

dev.off()


library(ggrepel)   # for nice non-overlapping labels

group_col <- "sctype_default"   # your metadata column
red <- "rna.umap"

# helper that adds labels to a FeaturePlot
featureplot_with_labels <- function(obj, feature, group_col = "sctype_default", reduction = "rna.umap") {
  # 1) Base plot: ADT overlay on RNA UMAP
  DefaultAssay(obj) <- "ADT"
  p <- FeaturePlot(
    obj,
    features = feature,
    reduction = reduction,
    min.cutoff = "q05", max.cutoff = "q95",
    cols = c("lightgrey", "blue")
  ) + ggtitle(paste0("Protein expression of ", feature, " on RNA UMAP"))
  
  # 2) Compute centroid per group on the same UMAP
  coords <- as.data.frame(Embeddings(obj, reduction = reduction)[, 1:2])
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  coords[[group_col]] <- obj@meta.data[[group_col]]
  
  cents <- coords %>%
    group_by(.data[[group_col]]) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
  
  # 3) Overlay labels
  p + ggrepel::geom_text_repel(
    data = transform(cents, label = cents[[group_col]]),
    aes(UMAP_1, UMAP_2, label = label),
    size = 4, fontface = "bold", segment.color = NA, box.padding = 0.25
  )
}

# Example for a single marker
featureplot_with_labels(seur_obj, "Hu.CD94")

# Example for a single marker
featureplot_with_labels(seur_obj, "ANXA2R") 


# multi-page PDF with labels for all ADT markers
pdf("ADT_markers_on_RNA_UMAP_with_labels2.pdf", width = 6, height = 5)
for (f in adt_markers) print(featureplot_with_labels(seur_obj, f, group_col, red))
dev.off()




for (f in adt_markers) {
  print(featureplot_with_labels(seur_obj, f, group_col, red))
}



DefaultAssay(seur_obj) <- "ADT"
adt_has_data <- ncol(GetAssayData(seur_obj, slot = "data")) > 0
if (!adt_has_data) {
  message("Normalizing ADT with CLR (data slot empty)...")
  seur_obj <- NormalizeData(seur_obj, normalization.method = "CLR")
  seur_obj <- ScaleData(seur_obj, verbose = FALSE)
}



saveRDS(seur_obj, file = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/harmony_22jan_annotated.rds")



# library(HGNChelper)
# 
# test <- c("CD56","CD16","CD94","NKp46","CD122","CD314","CD161","NKG2D")
# HGNChelper::checkGeneSymbols(test)
# 
# library(stringr)
# library(HGNChelper)
# 
# s <- "CD8,CD56,CD2,CD16,FCGR3A ,FCGR3B,CD94,CD3D,CD3E,CD3G,CD3Z,NKp46,CD11b,CD161,CD314,CD69,NKG7,CD122,NKG2D,GZMB,GZMA,GZMM,GNLY,COX6A2,ZMAT4,KIR2DL4"
# 
# tokens <- str_split(s, "\\s*,\\s*")[[1]]
# chk <- HGNChelper::checkGeneSymbols(tokens)
# 
# chk
# 

DefaultAssay(seur_obj) <- "RNA"

VlnPlot(seur_obj,
        features = c("TRAC","TRBC1","TRBC2","CD3D","CD3E"),
        group.by = "sctype_classification",
        pt.size = 0, ncol = 3)

# quantify % of NK cells expressing TCR transcripts
nk_cells <- WhichCells(seur_obj, expression = sctype_classification3 == "Natural killer  cells")
nk_meta <- FetchData(seur_obj, vars = c("TRAC","TRBC1","CD3D"), cells = nk_cells)

colMeans(nk_meta > 0)  # fraction of NK cells with detectable expression


nk_cells <- WhichCells(seur_obj, expression = sctype_classification == "Natural killer  cells")

mat_counts <- GetAssayData(seur_obj, slot = "counts", layer = "RNA")[c("TRAC","TRBC1","CD3D"), nk_cells, drop=FALSE]

# fraction of NK cells with >=1 UMI
colMeans(t(mat_counts) >= 1)




DefaultAssay(seur_obj) <- "RNA"

nk <- subset(seur_obj, subset = sctype_classification == "Natural killer  cells")

FeatureScatter(nk, feature1 = "TRBC1", feature2 = "TRAC")
FeatureScatter(nk, feature1 = "TRBC1", feature2 = "CD3D")
FeatureScatter(nk, feature1 = "TRAC",  feature2 = "CD3D")
FeatureScatter(nk, feature1 = "TRBC2", feature2 = "CD3E")


FeatureScatter(nk, feature1 = "KLRD1", feature2 = "TRBC2", slot = "counts")
FeatureScatter(nk, feature1 = "KLRD1", feature2 = "TRBC1", slot = "counts")

