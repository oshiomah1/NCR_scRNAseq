library(Seurat)

# Load your object
obj <- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/raw_merge_all_batches_harm_annotatedbatch1_2.rds")
DefaultAssay(obj) <- "RNA"


library(dplyr)
#downsample for testing 
# how many cells max per donor × cell type?
max_cells_per_donor_ct <- 50  # or 100, 500, etc.

meta <- obj@meta.data %>%
  mutate(cell = rownames(.))

# Split by donor × cell type, then sample within each group
keep_cells_tbl <- meta %>%
  group_split(NCR.ID, sctype_classification) %>%
  lapply(function(df) {
    n_keep <- min(nrow(df), max_cells_per_donor_ct)
    dplyr::slice_sample(df, n = n_keep)
  }) %>%
  bind_rows()

keep_cells <- keep_cells_tbl$cell

 

obj_thin <- subset(obj, cells = keep_cells)
obj_thin
rm(obj)

# Create a pseudobulk sample ID = donor × cell type
obj_thin$pb_id <- interaction(obj_thin$NCR.ID,
                              obj_thin$sctype_classification,
                              drop = TRUE)
# RNA pseudobulk
DefaultAssay(obj_thin) <- "RNA"
rna_pb <- AggregateExpression(
  obj_thin,
  assays   = "RNA",
  group.by = "pb_id",
  slot     = "counts",
  fun      = "sum"
)$RNA

# ADT pseudobulk
DefaultAssay(obj_thin) <- "ADT"
adt_pb <- AggregateExpression(
  obj_thin,
  assays   = "ADT",
  group.by = "pb_id",
  slot     = "counts",
  fun      = "sum"
)$ADT

dim(rna_pb)  # genes × pseudobulk samples
dim(adt_pb)  # proteins × pseudobulk samples
head(colnames(rna_pb))
head(colnames(adt_pb))  # should be identical


#Build the pseudobulk metadata (pb_meta)

pb_meta <- obj_thin@meta.data %>%
  mutate(pb_id = obj_thin$pb_id) %>%
  group_by(pb_id) %>%
  summarise(
    NCR.ID              = first(NCR.ID),
    sctype              = first(sctype_classification),
    validated_TB_status = first(validated_TB_status),
    Age                 = first(Age.),
    Gender              = first(Gender),
    group               = first(group),
    n_cells             = n(),
    .groups             = "drop"
  ) %>%
  # only keep those pb_ids that are actually in the count matrices
  filter(pb_id %in% colnames(rna_pb)) %>%
  arrange(pb_id)

# Make sure ordering matches
all(pb_meta$pb_id == colnames(rna_pb))

library(DESeq2)

# Choose a cell type
ct <- "ISG expressing immune cells"   # change as needed

# Subset metadata to that cell type
pb_meta_ct <- pb_meta %>%
  filter(sctype == ct)

# Require at least e.g. 2 cases and 2 controls
table(pb_meta_ct$validated_TB_status)

# Subset RNA counts to these pseudobulk samples
rna_ct_counts <- rna_pb[, pb_meta_ct$pb_id, drop = FALSE]

# Build DESeq2 object
dds_ct <- DESeqDataSetFromMatrix(
  countData = round(rna_ct_counts),
  colData   = as.data.frame(pb_meta_ct),
  design    = ~ validated_TB_status   # model: TB status effect per cell type
)

# Optional: filter lowly-expressed genes
keep <- rowSums(counts(dds_ct) >= 10) >= 3  # tweak thresholds as desired
dds_ct <- dds_ct[keep, ]

# Run DESeq2
dds_ct <- DESeq(dds_ct)

res_ct <- results(
  dds_ct,
  contrast = c("validated_TB_status", "2weekCase", "ctrl")
)

# Order by adjusted p-value
res_ct <- lfcShrink(dds_ct,
                    coef = "validated_TB_status_ctrl_vs_2weekCase",
                    type = "apeglm") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

head(res_ct)






######del below
# RNA
rna_pseudobulk <- AggregateExpression(
  obj_thin,
  assays  = "RNA",
  group.by = c("NCR.ID", "sctype_classification"),  
  slot   = "counts",                   # raw counts for DESeq2/edgeR-style pseudobulk
  fun    = "sum"
)$RNA

# ADT pseudobulk
 
adt_pseudobulk <- AggregateExpression(
  obj_thin,
  assays  = "ADT",
  group.by = c("NCR.ID", "sctype_classification"),  
  slot   = "counts",                   # raw counts for DESeq2/edgeR-style pseudobulk
  fun    = "sum"
)$ADT


dim(obj_thin)             # cells should be much smaller than 233k
dim(rna_pseudobulk)       # genes × pseudobulk samples
dim(adt_pseudobulk)       # proteins × same pseudobulk samples
head(colnames(rna_pseudobulk))
head(colnames(adt_pseudobulk))  # should match exactly


counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")

ncol(counts_mat) == ncol(obj)          # should be TRUE
all(colnames(counts_mat) == colnames(obj))  # should also be TRUE
tes = obj@meta.data


if ("ADT" %in% Assays(obj)) {
  obj[["ADT"]] <- JoinLayers(obj[["ADT"]])
}
##
 
library(SingleCellExperiment)

sce <- SingleCellExperiment(
  assays = list(counts = counts_mat),
  colData = obj@meta.data
)

library(scran)
#Pseudobulk aggregation (sample × celltype)
pb <- aggregateAcrossCells(
  sce,
  ids = DataFrame(
    sample   = sce$NCR.ID,
    celltype = sce$sctype_classification
  )
)

pb


library(edgeR)
library(dplyr)

pb_counts <- assay(pb, "counts")
pb_meta   <- as.data.frame(colData(pb))

 

# Quick sanity checks
table(pb_meta$celltype)[1:10]
table(pb_meta$validated_TB_status)        # should show Case vs Control (or similar)
table(pb_meta$celltype, pb_meta$validated_TB_status)[1:5, ]

#Run edgeR DGE per cell type (Case vs Control)

#Here’s a reusable loop that: subsets to one cell type,filters lowly-expressed genes,fits NB-GLM with TBstatus

#andstores full results table per cell type


 

 
pb_counts <- assay(pb, "counts")
pb_meta   <- as.data.frame(colData(pb))

# create a clean TBstatus variable from validated_TB_status
pb_meta$TBstatus <- factor(pb_meta$validated_TB_status)
pb_meta$TBstatus <- droplevels(pb_meta$TBstatus)
pb_meta$TBstatus <- relevel(pb_meta$TBstatus, ref = "ctrl")
pb_meta$batch <- factor(pb_meta$batch)
table(pb_meta$batch, pb_meta$TBstatus)
levels(pb_meta$batch)
 
levels(pb_meta$TBstatus)  # should be "ctrl" "2weekCase"

results_list <- list()

for (ct in sort(unique(pb_meta$celltype))) {
  message("Running DGE for cell type: ", ct)
  
  idx <- which(pb_meta$celltype == ct)
  
  # subset metadata and drop unused levels
  ct_meta <- droplevels(pb_meta[idx, , drop = FALSE])
  
  # skip if only one TBstatus present
  if (length(unique(ct_meta$TBstatus)) < 2) {
    message("  Skipping ", ct, " (only one TBstatus level present)")
    next
  }
  
  ct_counts <- pb_counts[, idx, drop = FALSE]
  
  # filter low-expressed genes
  keep <- rowSums(ct_counts > 5) >= 3
  ct_counts <- ct_counts[keep, , drop = FALSE]
  
  if (sum(keep) < 50) {
    message("  Skipping ", ct, " (too few genes after filtering)")
    next
  }
  
  dge <- DGEList(counts = ct_counts)
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~ TBstatus, data = ct_meta)
  message("  design columns: ", paste(colnames(design), collapse = ", "))
  # should be: (Intercept), TBstatus2weekCase
  
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  
  # use the 2nd column (TBstatus term) automatically
  lrt <- glmLRT(fit, coef = 2)
  
  res <- topTags(lrt, n = Inf)$table
  res$gene     <- rownames(res)
  res$celltype <- ct
  
  results_list[[ct]] <- res
}

length(results_list)


#######################
#edger batch adjusted
results_list_batchadj <- list()

for (ct in sort(unique(pb_meta$celltype))) {
  message("Running batch-adjusted DGE for cell type: ", ct)
  
  idx <- which(pb_meta$celltype == ct)
  ct_meta <- droplevels(pb_meta[idx, , drop = FALSE])
  
  # skip if only one TBstatus present
  if (length(unique(ct_meta$TBstatus)) < 2) {
    message("  Skipping ", ct, " (only one TBstatus level present)")
    next
  }
  
  ct_counts <- pb_counts[, idx, drop = FALSE]
  
  # filter low-expressed genes
  keep <- rowSums(ct_counts > 5) >= 3
  ct_counts <- ct_counts[keep, , drop = FALSE]
  if (sum(keep) < 50) {
    message("  Skipping ", ct, " (too few genes after filtering)")
    next
  }
  
  dge <- DGEList(counts = ct_counts)
  dge <- calcNormFactors(dge)
  
  #  design: batch + TBstatus
  design <- model.matrix(~ batch + TBstatus, data = ct_meta)
  message("  design columns: ", paste(colnames(design), collapse = ", "))
  # Typically: (Intercept), batchbatch2_paths, TBstatus2weekCase
  
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  
  # TBstatus should be the LAST column → use coef = ncol(design)
  lrt <- glmLRT(fit, coef = ncol(design))
  
  res <- topTags(lrt, n = Inf)$table
  res$gene     <- rownames(res)
  res$celltype <- ct
  
  results_list_batchadj[[ct]] <- res
}

length(results_list_batchadj)
#####################
###############################################################
#Combine and summarize results across cell types
###############################################################
 

all_res <- bind_rows(results_list)  # each element already had $celltype + $gene

# Mark significance
all_res <- all_res %>%
  mutate(
    is_sig = FDR < 0.05,
    direction = case_when(
      is_sig & logFC > 0 ~ "up_in_case",
      is_sig & logFC < 0 ~ "down_in_case",
      TRUE ~ "ns"
    )
  )

# How many DE genes per cell type?
deg_summary <- all_res %>%
  group_by(celltype) %>%
  summarise(
    n_sig = sum(is_sig),
    n_up  = sum(direction == "up_in_case"),
    n_down = sum(direction == "down_in_case"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig))

deg_summary





######## volcano plot 

ct <- "ISG expressing immune cells"

res_ct <- results_list[[ct]] %>%
  mutate(
    is_sig = FDR < 0.05,
    label = ifelse(is_sig & abs(logFC) > 1, gene, NA)
  )

head(res_ct, 20)

sig_ct <- res_ct %>%
  filter(FDR < 0.05)

nrow(sig_ct)
head(sig_ct, 20)

volcanoplot=ggplot(res_ct, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = is_sig), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = paste0(" DGE: ", ct, " (2weekCase vs ctrl)"),
    x = "log2 fold-change (case vs ctrl)",
    y = "-log10(p-value)",
    color = "FDR < 0.05"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/volcano_file.png",
  plot = volcanoplot,
  width = 6,
  height = 4.5,
  dpi = 300
)
######## heatmap of top genes

library(pheatmap)

# pick, say, top 30 genes by |logFC| in Classical Monocytes
top_genes_ct <- res_ct %>%
  arrange(FDR, desc(abs(logFC))) %>%
  filter(FDR < 0.05) %>%
  head(30) %>%
  pull(gene)

# extract those genes from the pseudobulk counts for that cell type
idx <- which(pb_meta$celltype == ct)
ct_counts <- assay(pb, "counts")[top_genes_ct, idx, drop = FALSE]

# log-transform for visualization
log_ct <- log2(ct_counts + 1)

ann <- pb_meta[idx, c("validated_TB_status","validated_TB_status", "NCR.ID")]
rownames(ann) <- colnames(log_ct)

heattmap=pheatmap(
  log_ct,
  
  show_rownames = TRUE,
  show_colnames = FALSE,
  clustering_distance_cols = "correlation",
  main = paste0("Top DE genes in ", ct, " (pseudobulk)")
)
pdf(file.path("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/process_with_wgs/", "heatmap_manual.pdf"), width = 7, height = 6)
grid::grid.newpage()
grid::grid.draw(heattmap$gtable)
dev.off()
 


# Join layers **on the assay**, NOT the whole object
#obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

#  Confirm layers are now collapsed
#Layers(obj[["RNA"]])
