

#read in meta data from suerat object and select relevant columns
 

meta <- seur_obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%   # ✅ add barcodes
  dplyr::select(
    cell,
    NCR.ID,
    validated_TB_status,
    Age.,
    Gender,
    sctype_classification_man, group
  ) %>%
  dplyr::filter(!is.na(NCR.ID))

# aggregate sub cell types into borader cell types
meta$celltype_simplified <- dplyr::case_when(
  meta$sctype_classification_man %in% c("Naive CD4+ T cells", "Memory CD4+ T cells") ~ "CD4+ T cells",
  meta$sctype_classification_man %in% c("Naive CD8+ T cells", "Effector CD8+ T cells", "CD8+ NKT-like cells") ~ "CD8+ T cells",
  meta$sctype_classification_man %in% c("Natural killer  cells") ~ "NK cells",
  meta$sctype_classification_man %in% c("Naive B cells", "Plasma B cells") ~ "B cells",
  meta$sctype_classification_man %in% c("Classical Monocytes", "Non-classical monocytes") ~ "Monocytes",
  meta$sctype_classification_man %in% c("Plasmacytoid Dendritic cells") ~ "Dendritic cells",
  TRUE ~ meta$sctype_classification_man
)

# pseudobulk function RNA
pseudobulk_rna <- function(
    obj,
    celltype,
    donor_col = "NCR.ID"
) {
  
  cells_use <- WhichCells(
    obj,
    expression = celltype_simplified == celltype
  )
  
  counts <- GetAssayData(
    obj,
    assay = "RNA",
    slot = "counts"
  )[, cells_use, drop = FALSE]
  
  donors <- obj@meta.data[cells_use, donor_col]
  
  pb_counts <- rowsum(
    t(counts),
    group = donors
  )
  
  t(pb_counts)
}


# EXample usage 
rna_pb_DC <- pseudobulk_rna(seur_obj, "Dendritic cells")


# ADT pseudobullk function
pseudobulk_adt <- function(
  obj,
  celltype,
  donor_col = "NCR.ID"
) {
  
  cells_use <- WhichCells(
    obj,
    expression = sctype_classification_man == celltype
  )
  
  counts <- GetAssayData(
    obj,
    assay = "ADT",
    slot = "counts"
  )[, cells_use, drop = FALSE]
  
  donors <- obj@meta.data[cells_use, donor_col]
  
  pb_counts <- rowsum(
    t(counts),
    group = donors
  )
  
  t(pb_counts)
}
