## =========================
## Simple line-by-line QC (single replicate)
## =========================

## ---- Libraries ----
library(Seurat)
library(ggplot2)
library(scales)
library(dplyr)


set.seed(411)
options(future.globals.maxSize = 3 * 1024^3)

## ---- Inputs you change here ----
replicate <- "TBSS01"   # <— change per run

## Base paths
base_cr <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq"
out_dir  <- file.path(base_cr, "results/seurat/1_QC_wgs/batch1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Files for this replicate
filtered_matrix_path <- file.path(base_cr, "results/cellranger/batch1", replicate,
                                  "per_sample_outs", replicate, "count",
                                  "sample_filtered_feature_bc_matrix.h5")
vireo_res_path       <- file.path(base_cr, "results/demultiplex/vireob1_with_WGS",
                                  replicate, "donor_ids.tsv")
tcr_prefix <- file.path(base_cr, "results/cellranger/batch1", replicate,
                        "per_sample_outs", replicate, "vdj_t/")
bcr_prefix <- file.path(base_cr, "results/cellranger/batch1", replicate,
                        "per_sample_outs", replicate, "vdj_b/")

## Outputs
pdf_path <- file.path(out_dir, paste0(replicate, "_basic_seurat.pdf"))
log_path <- file.path(out_dir, paste0(replicate, "_log.txt"))
rds_path <- file.path(out_dir, paste0(replicate, "qc_only_vdj.rds"))

## ---- Start logging & PDF ----
sink(file = log_path, split = TRUE)
on.exit({ try(dev.off(), silent=TRUE); try(sink(), silent=TRUE) }, add = TRUE)

cat("\n[INFO]", replicate, "- QC run started\n")
cat("[INFO] H5:", filtered_matrix_path, "\n")
cat("[INFO] Vireo:", vireo_res_path, "\n")
cat("[INFO] PDF:", pdf_path, "\n")
cat("[INFO] LOG:", log_path, "\n\n")

pdf(file = pdf_path, width = 8, height = 6)

## =========================
## STEP 1: Load 10x data
## =========================
testrun <- Read10X_h5(filtered_matrix_path)
cat("[OK] 10x H5 loaded\n")
str(testrun)

## Create Seurat object (RNA) and add ADT if present
if ("Gene Expression" %in% names(testrun)) {
  seurat_obj <- CreateSeuratObject(counts = testrun$`Gene Expression`, project = replicate)
  cat("[OK] Seurat object created (RNA)\n")
} else {
  stop("Gene Expression matrix not found in H5")
}
if ("Antibody Capture" %in% names(testrun)) {
  seurat_obj[["ADT"]] <- CreateAssayObject(counts = testrun$`Antibody Capture`)
  cat("[OK] ADT assay added\n")
} else {
  cat("[WARN] No ADT assay present\n")
}

## =========================
## STEP 2: Basic QC features
## =========================
seurat_obj[["percent.mt"]]   <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[LS]")
print(seurat_obj)

## =========================
## STEP 3: Platelet filter
## =========================
platelet_genes <- c("PF4","GNG11","PPBP")
present <- intersect(platelet_genes, rownames(seurat_obj))
cat("Present platelet genes:", paste(present, collapse=", "), "\n")
cat("Missing platelet genes:", paste(setdiff(platelet_genes, present), collapse=", "), "\n")
if (length(present) == 0) stop("No platelet genes found in RNA assay.")

seurat_obj$platelet_score <- Matrix::colSums(seurat_obj[["RNA"]]$counts[present, , drop=FALSE])
threshold <- 1
cat("Cells before platelet filter:", ncol(seurat_obj), "\n")
hist(seurat_obj$platelet_score, breaks = 5000, main="Platelet score", xlab="score",
     col="lightblue", border="black")
seurat_obj <- subset(seurat_obj, subset = platelet_score < threshold)
cat("Cells after platelet filter:", ncol(seurat_obj), "\n")
#backup = seurat_obj= backup
## =========================
## STEP 4: Add Vireo metadata
## =========================
vireo_res <- read.table(vireo_res_path, header = TRUE)
metadata = read.csv("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/metadata/batch1_metadata.csv", header = TRUE)
vireo_res_metadata = left_join(vireo_res,metadata, by = c("donor_id" = "Sample_ID"))
rownames(vireo_res_metadata) <- vireo_res_metadata$cell
## keep only cells present
keep_cells <- intersect(colnames(seurat_obj), rownames(vireo_res_metadata))
cat("Cells with Vireo metadata:", length(keep_cells), "\n")
seurat_obj <- subset(seurat_obj, cells = keep_cells)
seurat_obj <- AddMetaData(seurat_obj, vireo_res_metadata[keep_cells, c("donor_id","prob_max","NCR.ID" , "group", "Age." , "Gender" , "validated_TB_status", "Cell.Viability" ,"Live.Cell.Count" ,"other" )])

seurat_obj$singlet_source_vireo <- dplyr::case_when(
  seurat_obj$donor_id == "doublet"    ~ "DBL",
  seurat_obj$donor_id == "unassigned" ~ "unassigned",
  TRUE                                ~ "SNG"
)

## =========================
## STEP 5: Add VDJ clonotypes (T, then B)
## =========================
## helper: read and join clonotype annotations
read_clono <- function(prefix) {
  ann <- file.path(prefix, "filtered_contig_annotations.csv")
  clo <- file.path(prefix, "clonotypes.csv")
  if (!file.exists(ann) || !file.exists(clo)) return(NULL)
  tcr <- read.csv(ann)
  tcr <- tcr[!duplicated(tcr$barcode), c("barcode","raw_clonotype_id","v_gene")]
  names(tcr)[2] <- "clonotype_id"
  clono <- read.csv(clo)[, c("clonotype_id","cdr3s_aa")]
  merge(tcr, clono, by="clonotype_id", all.x=TRUE)
}

## TCR
t_dat <- read_clono(tcr_prefix)
if (!is.null(t_dat)) {
  rownames(t_dat) <- t_dat$barcode
  t_dat <- t_dat[, c("clonotype_id","v_gene","cdr3s_aa")]
  colnames(t_dat) <- paste0("t_", colnames(t_dat))
  seurat_obj <- AddMetaData(seurat_obj, t_dat[intersect(colnames(seurat_obj), rownames(t_dat)), , drop=FALSE])
  cat("[OK] TCR clonotypes added\n")
} else cat("[WARN] No TCR clonotype files found\n")

## BCR
b_dat <- read_clono(bcr_prefix)
if (!is.null(b_dat)) {
  rownames(b_dat) <- b_dat$barcode
  b_dat <- b_dat[, c("clonotype_id","v_gene","cdr3s_aa")]
  colnames(b_dat) <- paste0("b_", colnames(b_dat))
  seurat_obj <- AddMetaData(seurat_obj, b_dat[intersect(colnames(seurat_obj), rownames(b_dat)), , drop=FALSE])
  cat("[OK] BCR clonotypes added\n")
} else cat("[WARN] No BCR clonotype files found\n")

## =========================
## STEP 6: Keep singlets only
## =========================
seurat_obj <- subset(seurat_obj, subset = singlet_source_vireo == "SNG")
cat("Cells after keeping singlets:", ncol(seurat_obj), "\n")

## Snapshot counts/features
rna_cells <- ncol(seurat_obj[["RNA"]])
adt_cells <- if ("ADT" %in% names(seurat_obj@assays)) ncol(seurat_obj[["ADT"]]) else NA_integer_
cat("RNA cells:", rna_cells, " | ADT cells:", adt_cells, "\n")

## =========================
## STEP 7: MAD-based bounds & QC plots
## =========================
median_nFeature <- median(seurat_obj$nFeature_RNA)
mad_nFeature    <- mad(seurat_obj$nFeature_RNA)
lower_feat <- abs(median_nFeature - (3 * mad_nFeature))
upper_feat <- median_nFeature + (3 * mad_nFeature)

median_ncount_rna <- median(seurat_obj$nCount_RNA)
mad_ncount_rna    <- mad(seurat_obj$nCount_RNA)
upper_ncount_rna  <- median_ncount_rna + (3 * mad_ncount_rna)

if ("ADT" %in% names(seurat_obj@assays)) {
  median_nFeature_adt <- median(seurat_obj$nFeature_ADT)
  mad_nFeature_adt    <- mad(seurat_obj$nFeature_ADT)
  upper_feat_adt      <- median_nFeature_adt + 3 * mad_nFeature_adt
  
  median_ncount_adt <- median(seurat_obj$nCount_ADT)
  mad_ncount_adt    <- mad(seurat_obj$nCount_ADT)
  upper_ncount_adt  <- median_ncount_adt + 3 * mad_ncount_adt
} else {
  upper_feat_adt <- Inf; upper_ncount_adt <- Inf
}

cat(sprintf("nFeature_RNA median=%g MAD=%g → [%g, %g]\n",
            median_nFeature, mad_nFeature, lower_feat, upper_feat))
cat(sprintf("nCount_RNA  median=%g MAD=%g → upper=%g\n",
            median_ncount_rna, mad_ncount_rna, upper_ncount_rna))
cat(sprintf("nFeature_ADT upper≈%g; nCount_ADT upper≈%g\n", upper_feat_adt, upper_ncount_adt))

print(VlnPlot(seurat_obj, features = "nFeature_RNA", group.by = "donor_id") + guides(fill="none") +
        geom_hline(yintercept = upper_feat))
print(VlnPlot(seurat_obj, features = "nCount_RNA", group.by = "donor_id") + guides(fill="none") +
        geom_hline(yintercept = upper_ncount_rna))
print(FeatureScatter(seurat_obj, feature1="nCount_RNA", feature2="percent.mt") +
        geom_hline(yintercept=15) + geom_density_2d())
print(FeatureScatter(seurat_obj, feature1="nFeature_RNA", feature2="percent.mt") +
        geom_vline(xintercept=lower_feat) + geom_vline(xintercept=250) +
        geom_vline(xintercept=upper_feat) + geom_density_2d())
print(VlnPlot(seurat_obj, features="percent.mt", group.by="donor_id") + guides(fill="none"))

if ("ADT" %in% names(seurat_obj@assays)) {
  print(VlnPlot(seurat_obj, features="nCount_ADT", group.by="donor_id") + guides(fill="none") +
          geom_hline(yintercept=upper_ncount_adt))
  print(VlnPlot(seurat_obj, features="nFeature_ADT", group.by="donor_id") + guides(fill="none") +
          geom_hline(yintercept=upper_feat_adt))
  print(FeatureScatter(seurat_obj, feature1="nCount_ADT", feature2="nFeature_ADT"))
}
print(VlnPlot(seurat_obj, features="percent.ribo", group.by="donor_id") + guides(fill="none"))

## =========================
## STEP 8: Final filtering
## =========================
seurat_obj <- subset(
  seurat_obj,
  subset =
    nFeature_RNA > lower_feat &
    nFeature_RNA < upper_feat &
    percent.mt < 15 &
    nCount_ADT < upper_ncount_adt
)
cat("Cells after final filter:", ncol(seurat_obj), "\n")

print(FeatureScatter(seurat_obj, feature1="nCount_RNA", feature2="nFeature_RNA"))

## =========================
## STEP 9: Clonotype summaries
## =========================
cat("Unique TCR CDR3s:", length(unique(seurat_obj$t_cdr3s_aa)), "\n")
cat("Unique BCR CDR3s:", length(unique(seurat_obj$b_cdr3s_aa)), "\n")
print(table(!is.na(seurat_obj$t_clonotype_id), !is.na(seurat_obj$b_clonotype_id)))

tcr_freq <- sort(table(seurat_obj$t_cdr3s_aa), decreasing=TRUE)
bcr_freq <- sort(table(seurat_obj$b_cdr3s_aa), decreasing=TRUE)
bcr_freq_df <- as.data.frame(bcr_freq)

top_10_t <- head(tcr_freq[tcr_freq >= 10], 10)
top_10_b <- head(bcr_freq[bcr_freq >= 10], 10)
print(top_10_t); print(top_10_b)

if (length(top_10_t) > 0) {
  tdf <- data.frame(clonotype = names(top_10_t), frequency = as.vector(top_10_t))
  print(
    ggplot(tdf, aes(x = reorder(clonotype, -frequency), y = frequency)) +
      geom_bar(stat="identity", fill="#B00B69") + coord_flip() +
      labs(title="TCR Clonotype Frequency", x="Clonotype", y="Frequency") +
      theme_minimal() +
      theme(axis.text.x=element_text(size=5, face="bold"),
            axis.text.y=element_text(size=5, face="bold"),
            plot.title=element_text(size=5, face="bold"))
  )
}
if (length(top_10_b) > 0) {
  bdf <- data.frame(clonotype = names(top_10_b), frequency = as.vector(top_10_b))
  print(
    ggplot(bdf, aes(x = reorder(clonotype, -frequency), y = frequency)) +
      geom_bar(stat="identity", fill="#69B00B") + coord_flip() +
      labs(title="BCR Clonotype Frequency", x="Clonotype", y="Frequency") +
      theme_minimal() +
      theme(axis.text.x=element_text(size=5, face="bold"),
            axis.text.y=element_text(size=5, face="bold"),
            plot.title=element_text(size=5, face="bold"))
  )
}

## =========================
## STEP 10: Save object & finish
## =========================
#saveRDS(seurat_obj, file = rds_path)
cat("\n[DONE] Saved:", rds_path, "\n")
cat("[DONE] PDF:", pdf_path, "\n")
#dev.off()
