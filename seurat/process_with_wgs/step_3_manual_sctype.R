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
seur_obj <- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/raw_merge_all_batches_harm_annotated_all_res6_pca15.rds") #reg processing + wnn

print("read file")
head(seur_obj)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# DB file
#db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx";
db_ <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/ScTypeDB_immune_only.xlsx"
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

DimPlot(seur_obj, reduction = "rna.umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification_man')        

saveRDS(seur_obj, file = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/raw_merge_all_batches_harm_sc_annotated_all_res6_pca15.rds")

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

mygraph <- graph_from_data_frame(edges, vertices=nodes)
# Make the graph
gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

DimPlot(seur_obj, reduction = "rna.umap", label = TRUE, repel = TRUE, cols = ccolss)+ gggr
DimPlot(seur_obj, reduction = "rna.umap", label = TRUE, repel = TRUE)+ gggr
library(pheatmap)

# Build matrix: rows = cell types, cols = clusters
mat <- cL_resutls %>%
  dplyr::select(cluster, type, scores) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = scores, values_fill = 0) %>%
  tibble::column_to_rownames("type") %>%
  as.matrix()

pheatmap(
  mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  fontsize_col = 9
)

top4 <- cL_resutls %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(scores, n = 4, with_ties = FALSE) %>%
  dplyr::ungroup()
library(ggplot2)


pdf("ScType_top4_per_cluster_pcs15_new.pdf", width = 7, height = 5)
 

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

