# Load required libraries
library(Seurat)
library(ggplot2)
library(scales)
library(dplyr)
#library(DoubletFinder)
library(glmGamPoi)
library(future)
library(furrr)

# Set the seed for reproducibility
set.seed(411)
print("m")
options(future.globals.maxSize = 3 * 1024^3)

# Replicates to process
 

replicates <- c("TBSS01") #here

# Output PDF path (will be updated for each replicate)
output_pdf_base <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_test/batch1/" #here
# Output log file path (will be updated for each replicate)
log_file_base <- "/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_test/batch1/"   #here

# Set up parallelization plan
plan(multisession, workers = 1)  # Adjust the number of workers based on your system

# Function to process each replicate
process_replicate <- function(replicate) {
  # Set the output PDF file path
  print(paste("Processing replicate:", replicate))
  
  # Set paths dynamically for each replicate
  output_pdf <- paste0(output_pdf_base, replicate, "_basic_seurat.pdf")
  log_file <-  paste0(log_file_base, replicate, "_log.txt")
  print(log_file)
  print(output_pdf)
  # Redirect console output to the log file
  sink(file = log_file, split = TRUE)
  print("Logging started...")
  filtered_matrix_path <- paste0("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/cellranger/batch1/", replicate, "/per_sample_outs/", replicate, "/count/sample_filtered_feature_bc_matrix.h5") #here
  vireo_res_path <- paste0("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/demultiplex/vireob1_with_WGS/", replicate, "/donor_ids.tsv") #here
  
  # Open the PDF device to save plots
  pdf(file = output_pdf, width = 8, height = 6)
  
  
  # Read data from HDF5
  testrun <- Read10X_h5(filtered_matrix_path)
  print("File is loaded")
  
  # Check the structure of the data
  str(testrun)
  
  # Create Seurat object for the gene expression data
  if ("Gene Expression" %in% names(testrun)) {
    seurat_obj <- CreateSeuratObject(counts = testrun$`Gene Expression`, project = replicate) #here
    print(paste("Seurat object created for Gene Expression data -", replicate))
  }
  
  # Add antibody capture data as a new assay, if present
  if ("Antibody Capture" %in% names(testrun)) {
    seurat_obj[["ADT"]] <- CreateAssayObject(counts = testrun$`Antibody Capture`)
    print("Antibody Capture data added as ADT assay")
  }
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  print("Calculated mitochondrial percentage")
  
  # Check the Seurat object
  print(seurat_obj)
  
  
  #grep("^RP[LS]",rownames(data),value = TRUE)
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[LS]")
  
  
  # List of platelet-specific genes
  platelet_genes <- c("PF4", "GNG11", "PPBP")
  
  # Check if the genes exist in the dataset
  present_genes <- intersect(platelet_genes, rownames(seurat_obj))
  
  # Display which genes are present and missing
  cat("Present genes:", present_genes, "\n")
  cat("Missing genes:", setdiff(platelet_genes, present_genes), "\n")
  
  # If no platelet-specific genes are found, stop the process
  if (length(present_genes) == 0) {
    stop("None of the platelet-specific genes are present in the Seurat object.")
  }
  
  # Calculate the total expression of the present platelet-specific genes
  seurat_obj$platelet_score <- Matrix::colSums(seurat_obj[["RNA"]]$counts[present_genes, ])
  
  # Set a threshold to define platelet-like cells
  threshold <- 1  # Adjust this based on your data
  print(paste("Number of cells original f1:", ncol(seurat_obj)))
  # Plot histogram of platelet_score
  hist(seurat_obj@meta.data$platelet_score,
       breaks = 5000,  # Adjust the number of bins if needed
       main = "Histogram of Platelet Scores",
       xlab = "Platelet Score",
       col = "lightblue", 
       border = "black")
  
  seurat_obj <- subset(seurat_obj, subset = platelet_score < threshold)
  print(paste("Number of cells after platelet  filter(2):", ncol(seurat_obj)))
  
  # Check the added metadata
  
  str(seurat_obj)
  
  # Load Vireo results
  vireo_res <- read.table(vireo_res_path, header = TRUE)
  
  # Add the barcode as rownames to the vireo results
  rownames(vireo_res) <- vireo_res$cell
  length(rownames(vireo_res))
  # Extract the columns you want to add to Seurat metadata (e.g., SNG.1ST and DBL.1ST)
  seurat_obj <- AddMetaData(seurat_obj, vireo_res[, c("cell", "donor_id", "prob_max")]) 
  
  #make a new variable to record the three classes
  seurat_obj$singlet_source_vireo <- case_when(
    seurat_obj$donor_id == "doublet" ~ "DBL",
    seurat_obj$donor_id == "unassigned" ~ "unassigned",
    TRUE ~ "SNG"
  )
  
  add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
    tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))
    # ADD VGENE COLUNE L8R!!!!!
    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    tcr <- tcr[!duplicated(tcr$barcode), ]
    
    # Only keep the barcode and clonotype columns. 
    # We'll get additional clonotype info from the clonotype table.
    tcr <- tcr[,c("barcode", "raw_clonotype_id", "v_gene")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
    
    # Clonotype-centric info.
    clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))
    
    # Slap the AA sequences onto our original table by clonotype_id.
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
    names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
    print(head(tcr))
    # Reorder so barcodes are first column and set them as rownames.
    tcr <- tcr[, c(2,1,3,4)]
    print(head(tcr))
    rownames(tcr) <- tcr[,1]
    tcr[,1] <- NULL
    colnames(tcr) <- paste(type, colnames(tcr), sep="_")
    print(head(tcr))
    # Add to the Seurat object's metadata.
    clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
    return(clono_seurat)
  }
  
  seurat_obj <- add_clonotype(paste0("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/cellranger/batch1/" , replicate, "/per_sample_outs/", replicate, "/vdj_t/"), seurat_obj, "t")# here
  seurat_obj <- add_clonotype(paste0("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/cellranger/batch1/" , replicate, "/per_sample_outs/", replicate, "/vdj_b/"), seurat_obj, "b")# here
  
  head(seurat_obj@meta.data)
  
  # Filter out doublets
  
  seurat_obj <- subset(seurat_obj, subset = singlet_source_vireo == "SNG")
  
  #nrow(seurat_obj@meta.data) #no of cells after filtER for only singlets
  print(paste("Number of cells after doublets  filter(3):", ncol(seurat_obj)))
  
  # Check the ADT assay
  seurat_obj[["ADT"]]
  # List all ADT features
  length(rownames(seurat_obj[["ADT"]]))
  
  # Check the RNA assay
  seurat_obj[["RNA"]]
  # List all the RNA features
  length(rownames(seurat_obj[["RNA"]]))
  
  print(paste("Number of cells before filtering:", ncol(seurat_obj)))
  rna_cells <- ncol(seurat_obj[["RNA"]])
  adt_cells <- ncol(seurat_obj[["ADT"]])
  
  print(paste("Number of RNA cells:", rna_cells))
  print(paste("Number of ADT cells:", adt_cells))
  
  if (rna_cells == adt_cells) {
    print("The RNA and ADT assays have the same number of cells.")
  } else {
    print("The RNA and ADT assays have different numbers of cells.")
  }
  
  # Inspect RNA assay counts matrix
  rna_counts <- seurat_obj[["RNA"]]@layers[["counts"]]
  dim(rna_counts)  # Dimensions of the RNA count matrix
  
  
  # Inspect ADT assay counts matrix
  adt_counts <- seurat_obj[["ADT"]]@counts
  dim(adt_counts)  # Dimensions of the ADT count matrix
  
  # Calculate the median and MAD for nFeature_RNA
  median_nFeature <- median(seurat_obj$nFeature_RNA)
  mad_nFeature <- mad(seurat_obj$nFeature_RNA)
  print(paste("Median nFeature_RNA:", median_nFeature))
  print(paste("MAD nFeature_RNA:", mad_nFeature))
  
  # Define upper and lower bounds using MAD for nFeature_RNA
  lower_bound_rna_nfeature <- abs(median_nFeature - (3 * mad_nFeature))
  upper_bound_rna_nfeature <- median_nFeature + (3 * mad_nFeature)
  print(paste("Lower bound for nFeature_RNA (3*MAD):", lower_bound_rna_nfeature))
  print(paste("Upper bound for nFeature_RNA (3*MAD):", upper_bound_rna_nfeature))
  
  # Calculate the median and MAD for nFeature_ADT
  median_nFeature_adt <- median(seurat_obj$nFeature_ADT)
  mad_nFeature_adt <- mad(seurat_obj$nFeature_ADT)
  print(paste("Median nFeature_ADT:", median_nFeature_adt))
  print(paste("MAD nFeature_ADT:", mad_nFeature_adt))
  
  # Define upper and lower bounds using MAD for nFeature_ADT
  lower_bound_nFeature_adt <- abs(median_nFeature_adt - (3 * mad_nFeature_adt))
  upper_bound_nFeature_adt <- median_nFeature_adt + (3 * mad_nFeature_adt)
  print(paste("Lower bound for nFeature_ADT (3*MAD):", lower_bound_nFeature_adt))
  print(paste("Upper bound for nFeature_ADT (3*MAD):", upper_bound_nFeature_adt))
  
  # Calculate the median and MAD for nCount_ADT
  median_ncount_adt <- median(seurat_obj$nCount_ADT)
  mad_ncount_adt <- mad(seurat_obj$nCount_ADT)
  print(paste("Median nCount_ADT:", median_ncount_adt))
  print(paste("MAD nCount_ADT:", mad_ncount_adt))
  
  # Define upper and lower bounds using MAD for nCount_ADT
  lower_bound_ncount_adt <- abs(median_ncount_adt - (3 * mad_ncount_adt))
  upper_bound_ncount_adt <- median_ncount_adt + (3 * mad_ncount_adt)
  print(paste("Lower bound for nCount_ADT (3*MAD):", lower_bound_ncount_adt))
  print(paste("Upper bound for nCount_ADT (3*MAD):", upper_bound_ncount_adt))
  
  # Calculate the median and MAD for nCount_RNA
  median_ncount_rna <- median(seurat_obj$nCount_RNA)
  mad_ncount_rna <- mad(seurat_obj$nCount_RNA)
  print(paste("Median nCount_RNA:", median_ncount_rna))
  print(paste("MAD nCount_RNA:", mad_ncount_rna))
  
  # Define upper and lower bounds using MAD for nCount_RNA
  lower_bound_ncount_rna <- abs(median_ncount_rna - (3 * mad_ncount_rna))
  upper_bound_ncount_rna <- median_ncount_rna + (3 * mad_ncount_rna)
  print(paste("Lower bound for nCount_RNA (3*MAD):", lower_bound_ncount_rna))
  print(paste("Upper bound for nCount_RNA (3*MAD):", upper_bound_ncount_rna))
  
  #violin plots
  
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA"), group.by = "donor_id") + guides(fill = "none") + geom_hline(yintercept =upper_bound_rna_nfeature ))
  print(VlnPlot(seurat_obj, features = c("nCount_RNA"), group.by = "donor_id") + guides(fill = "none") + geom_hline(yintercept=upper_bound_ncount_rna ))
  print(FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_hline(yintercept = 15) + geom_density_2d())
  print(FeatureScatter(seurat_obj, feature1 = "nFeature_RNA", feature2 = "percent.mt") + geom_vline(xintercept = lower_bound_rna_nfeature) + geom_vline(xintercept = 250) +geom_vline(xintercept= upper_bound_rna_nfeature)+ geom_density_2d())
  print(VlnPlot(seurat_obj, features = c("percent.mt"), group.by = "donor_id") + guides(fill = "none"))
  print(VlnPlot(seurat_obj, features = c("nCount_ADT"), group.by = "donor_id") + guides(fill = "none") + geom_hline(yintercept=upper_bound_ncount_adt))
  print(VlnPlot(seurat_obj, features = c("nFeature_ADT"), group.by = "donor_id") + guides(fill = "none") + geom_hline(yintercept = upper_bound_ncount_adt))
  print(FeatureScatter(seurat_obj, feature1 = "nCount_ADT", feature2 = "nFeature_ADT"))
  print(VlnPlot(seurat_obj, features = c("percent.ribo"), group.by = "donor_id") + guides(fill = "none"))
  
  # Filter the cells based on feature counts and mitochondrial content
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > lower_bound_rna_nfeature & nFeature_RNA < upper_bound_rna_nfeature & percent.mt < 15 &  nCount_ADT < upper_bound_ncount_adt)
  print("Filtered the cells based on feature counts and mitochondrial content")
  print(paste("Number of cells after final  filter(4):", ncol(seurat_obj)))
  
  # Check the RNA assay after filtering
  seurat_obj[["RNA"]]
  # List all the RNA features
  length(rownames(seurat_obj[["RNA"]]))
  seurat_obj[["ADT"]]
  # List all the RNA features
  length(rownames(seurat_obj[["ADT"]]))
  
  #plot features post rna ASSAY  filtering
  print(FeatureScatter(seurat_obj, feature1 ="nCount_RNA", feature2 = "nFeature_RNA" ) )
  #+ geom_smooth(method ="lm")
  print("plotted features post filtering")
  print(paste("Number of unique TCR AFTER BASIC filtering:", length(unique(seurat_obj$t_cdr3s_aa))))
  print(paste("Number of unique BCR AFTER BASIC filtering:",length(unique(seurat_obj$b_cdr3s_aa))))
  #seurat_obj@meta.data
  table(!is.na(seurat_obj$t_clonotype_id),!is.na(seurat_obj$b_clonotype_id))
  #end
  # TCR clonotype frequency
  tcr_freq <- table(seurat_obj$t_cdr3s_aa)
  #tcr_freq <- table(seurat_obj$t_v_gene)
  # Display the top 10 clonotypes by frequency
  head(sort(tcr_freq, decreasing = TRUE), 10)
  
  # BCR clonotype frequency
  bcr_freq <- table(seurat_obj$b_cdr3s_aa)
  #bcr_freq <- table(seurat_obj$b_v_gene)
  
  # Display the top 10 clonotypes by frequency
  head(sort(bcr_freq, decreasing = TRUE), 10)
  
  # Identify expanded TCR clonotypes (appears in 10 or more cells)
  expanded_tcr <- tcr_freq[tcr_freq >= 10]
  # Get the top 10 expanded TCR clonotypes by frequency
  top_10_expanded_tcr <- head(sort(expanded_tcr, decreasing = TRUE), 10)
  # Display expanded TCR clonotypes
  length(expanded_tcr)
  length(top_10_expanded_tcr)
  print(top_10_expanded_tcr)
  # Identify expanded BCR clonotypes (appears in 10 or more cells)
  expanded_bcr <- bcr_freq[bcr_freq >= 10]
  top_10_expanded_bcr <- head(sort(expanded_bcr, decreasing = TRUE), 10)
  # Display expanded TCR clonotypes
  length(expanded_bcr)
  length(top_10_expanded_bcr)
  print(top_10_expanded_bcr)
  # Bar plot of TCR clonotype frequencies
  
  tcr_freq_df <- data.frame(clonotype = names(top_10_expanded_tcr), frequency = as.vector(top_10_expanded_tcr))
  
  tcrplot<- ggplot(tcr_freq_df, aes(x = reorder(clonotype, -frequency), y = frequency)) +
    geom_bar(stat = "identity", fill = "#B00B69") +
    coord_flip() + 
    labs(title = "TCR Clonotype Frequency", x = "Clonotype", y = "Frequency") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 5, face = "bold"), # Make x-axis labels bold and size 12
      axis.text.y = element_text(size = 5, face = "bold"),  # Make y-axis labels bold and size 12
      plot.title = element_text(size = 5, face = "bold")
    )
  # Bar plot of BCR clonotype frequencies
  
  bcr_freq_df <- data.frame(clonotype = names(top_10_expanded_bcr), frequency = as.vector(top_10_expanded_bcr))
  
  bcrplot<-ggplot(bcr_freq_df, aes(x = reorder(clonotype, -frequency), y = frequency)) +
    geom_bar(stat = "identity", fill = "#69B00B") +
    coord_flip() + 
    labs(title = "BCR Clonotype Frequency", x = "Clonotype", y = "Frequency") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 5, face = "bold"), # Make x-axis labels bold and size 12
      axis.text.y = element_text(size = 5, face = "bold"),
      plot.title = element_text(size = 5, face = "bold")
    )
  
  print(tcrplot)
  print(bcrplot)
  
  # Save the Seurat object as an RDS file #here
  saveRDS(seurat_obj, file = paste0("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/1_QC_test/batch1/", replicate, "qc_only_vdj.rds"))#here 
  print(paste("Processing complete for replicate:", replicate))
  
  dev.off()
  # Stop logging and restore console output
  sink()
  
}

# Run the processing for all replicates in parallel
future_map(replicates, process_replicate,.options = furrr_options(seed = TRUE))



pdf("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/seurat/test.pdf", width=8, height=6)
plot(1:10)
dev.off()

