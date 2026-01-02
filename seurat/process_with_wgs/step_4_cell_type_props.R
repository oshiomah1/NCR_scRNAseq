library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

# Load your seur_object
#obj <- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/3_harmonize_batches_wgs/raw_merge_all_batches_harm_annotatedbatch1_2.rds")
seur_obj <- readRDS("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/seurat/raw_merge_all_batches_harm_sc_annotated_all_res6_pca15.rds") #reg processing + wnn

# #RE ANNOTATE AFTER CHANGING TO CD 16
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R")
# seur_obj <- run_sctype(
#   seur_obj , assay="RNA", scaled=TRUE,
#   known_tissue_type="Immune system",
#   custom_marker_file="/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/ScTypeDB_shortCD16_TRAC.xlsx",
#   name="sctype_classification3"
# )



DefaultAssay(seur_obj) <- "RNA"
head(seur_obj@meta.data[, c("donor_id", "validated_TB_status", "sctype_classification_man")])
table(seur_obj$validated_TB_status)
table(seur_obj$donor_id)
table(seur_obj$sctype_classification_man)


library(dplyr)

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

 
# filter to an cell type as a test
df_plot <- celltype_props %>%
  filter(sctype_classification_man == "Natural killer  cells")

ggplot(
  df_plot,
  aes(
    x = validated_case_status,
    #x = validated_TB_status,
    y = prop,
    color = days_antibiotics
  )
) +
  geom_quasirandom(width = 0.2, size = 2.2, alpha = 0.85) +
  scale_color_viridis_c(
    option = "plasma",
    na.value = "grey70"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = NULL,
    y = "Proportion of NK cells",
    color = "Days on antibiotics",
    title = "NK cell proportions by case status"
  )

#plot testt cell type bbeswarm colorcoded by gender
ggplot(
  df_plot,
  aes(
    x = validated_TB_status,
    y = prop,
    color = Gender
  )
) +
  geom_quasirandom(width = 0.2, size = 2.2, alpha = 0.85) +
  scale_color_viridis_d(
    option = "plasma",
    na.value = "grey70"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = NULL,
    y = "Proportion of NK cells",
    color = "NCR_Gender",
    title = "NK cell proportions by case status"
  )

celltypes <- sort(unique(celltype_props$sctype_classification_man))

#function to automate beeswarm 
plot_celltype_beeswarm <- function(
    data,
    celltype,
    color_var,
    color_label,
    plot_title_suffix = "",
    continuous_palette = c("#2c7bb6", "#feffbf", "#d7191c"),
    discrete_palette = "Set2"
) {
  
  df_plot <- data %>%
    dplyr::filter(sctype_classification_man == celltype)
  
  is_continuous <- is.numeric(df_plot[[color_var]]) || is.integer(df_plot[[color_var]])
  
  color_scale <- if (is_continuous) {
    scale_color_gradientn(
      colors = continuous_palette,
      na.value = "grey70"
    )
  } else {
    scale_color_brewer(
      palette = discrete_palette,
      na.value = "grey70"
    )
  }
  
  ggplot(
    df_plot,
    aes(
      x = validated_TB_status,
      y = prop
    )
  ) +
    geom_boxplot(
      width = 0.35,
      outlier.shape = NA,
      fill = "grey85",
      alpha = 0.35,
      color = "black"
    ) +
    ggbeeswarm::geom_quasirandom(
      aes(color = .data[[color_var]]),
      width = 0.2,
      size = 2.2,
      alpha = 0.85
    ) +
    color_scale +
    theme_minimal(base_size = 12) +
    labs(
      x = NULL,
      y = paste("Proportion of", celltype),
      color = color_label,
      title = paste(celltype, plot_title_suffix)
    )
}

####################################
#########test run ###########
####################################
plot_celltype_beeswarm(
  data = df_plot,
  celltype = "Natural killer  cells",
  color_var = "Age.",
  color_label = "Age",
  continuous_palette = c("navy", "white", "firebrick")
)

plot_celltype_beeswarm(
  data = df_plot,
  celltype = "Natural killer  cells",
  color_var = "Gender",
  color_label = "Sex"
)
####################################

####################################
###GeNerate Age Beeswarm ##########
####################################

pdf("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/cell_type_proportions/Beeswarm_all_celltypes_age.pdf", width = 6, height = 4)

for (ct in celltypes) {
  p <- plot_celltype_beeswarm(
    data = celltype_props,
    celltype = ct,
    color_var = "Age.",
    color_label = "Age",
    continuous_palette = c("navy", "white", "firebrick")
  )
  print(p)
}

dev.off()


####################################
###GeNerate abxBeeswarm ##########
####################################

pdf("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/cell_type_proportions/Beeswarm_all_celltypes_days_antibiotics.pdf", width = 6, height = 4)

for (ct in celltypes) {
  p <- plot_celltype_beeswarm(
    data = celltype_props,
    celltype = ct,
    color_var = "days_antibiotics",
    color_label = "Days on antibiotics",
    continuous_palette = c("navy", "white", "firebrick")
  )
  print(p)
}

dev.off()

####################################
###GeNerate gender Beeswarm ##########
####################################
pdf("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/cell_type_proportions/Beeswarm_all_celltypes_gender.pdf", width = 6, height = 4)

for (ct in celltypes) {
  p <- plot_celltype_beeswarm(
    data = celltype_props,
    celltype = ct,
    color_var = "Gender",
    color_label = "Gender",
    
  )
  print(p)
}

dev.off()


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


############
#poxsibly delte below
#########
#   sctype vs TB status — bar charts (counts + proportions)

meta <- seur_obj@meta.data %>%
  transmute(
    sctype = sctype_classification,
    TBstatus = validated_TB_status
  ) %>%
  filter(!is.na(sctype), !is.na(TBstatus))

meta2 <- seur_obj@meta.data %>%
  transmute(donor_id, sctype = sctype_classification3,
            TBstatus = validated_TB_status) %>%
  filter(!is.na(donor_id), !is.na(sctype), !is.na(TBstatus)) %>%
  mutate(TBstatus = recode(TBstatus, "2weekCase"="Case", "ctrl"="Control"))



 

 

ct <- meta %>%
  dplyr::count(
    .data$TBstatus,   # group by the TBstatus column in `meta`
    .data$sctype,     # group by the sctype column in `meta`
    name = "cells"    # name the count column "cells"
  )

ct_donor <- meta2 %>%
  count(donor_id, TBstatus, sctype, name="cells") %>%
  group_by(donor_id, TBstatus) %>%
  mutate(prop = cells / sum(cells)) %>%
  ungroup()



 
# Order cell types by overall size
sctype_order <- ct %>% group_by(sctype) %>% summarise(tot = sum(cells)) %>%
  arrange(desc(tot)) %>% pull(sctype)
ct$sctype <- factor(ct$sctype, levels = sctype_order)

sctype_order2 <- ct_donor %>% group_by(sctype) %>% summarise(tot = sum(cells)) %>%
  arrange(desc(tot)) %>% pull(sctype)
ct_donor$sctype <- factor(ct_donor$sctype, levels = sctype_order)

# (a) Stacked counts per TB status
# print(
#   ggplot(ct, aes(x = TBstatus, y = cells, fill = sctype)) +
#     geom_col() +
#     labs(title = "Cell-type composition by validated_TB_status",
#          x = "validated_TB_status", y = "Cells", fill = "ScType") +
#     theme_minimal() +
#     theme(legend.position = "right")
# )
# 
# (b) Proportions within each TBstatus
ct_prop <- ct %>%
  group_by(TBstatus) %>%
  mutate(prop = cells / sum(cells)) %>%
  ungroup()

# print(
#   ggplot(ct_prop, aes(x = TBstatus, y = prop, fill = sctype)) +
#     geom_col() +
#     scale_y_continuous(labels = scales::percent_format()) +
#     labs(title = " ",
#          x = "TB status", y = " Cell Type Proportion", fill = "ScType") +
#     theme_minimal()
# )

# ---- Cell-ready stacked proportions: recode, ascending legend, nice styling ----
 
library(forcats)
 
library(scales)
library(grid)   # for unit()

# 1) Start from your ct_prop (cols: TBstatus, sctype, prop)

# Recode x-axis labels robustly
ct_prop2 <- ct_prop %>%
  mutate(
    TBstatus = factor(
      TBstatus,
      levels = c("2weekCase", "ctrl"),
      labels = c("Case", "Control")
    )
  )

# Order legend (fill) by ASCENDING overall abundance
order_fill <- ct_prop2 %>%
  group_by(sctype) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_prop)) %>%
  pull(sctype) %>%
  as.character()    


ct_prop2 <- ct_prop2 %>%
  mutate(sctype = fct_relevel(sctype, order_fill))

# Build a palette long enough for all levels (colorblind-friendly-ish default)
pal <- setNames(hue_pal()(nlevels(ct_prop2$sctype)), levels(ct_prop2$sctype))

# 2) Plot
                
okabe_ito <- c(
  "#000000", # black (anchor / outline)
  "#E69F00", # orange (Okabe–Ito)
  "#56B4E9", # light blue (keep one)
  "#009E73", # green
  "#F0E442", # yellow
  "#D55E00", # red-orange
  "#CC79A7", # magenta
  "#999999", # gray
  "#800080", # deep purple
  "#8B4513", # saddle brown
  "#20B2AA", # light sea green
  "#FFD700", # golden yellow
  "#DC143C", # crimson
  "#708090", # slate gray
  "#8A2BE2", # blue violet
  "#A52A2A"  # brick brown
)
p <- ggplot(ct_prop2, aes(x = TBstatus, y = prop, fill = sctype)) +
  geom_col(position = position_stack(reverse = TRUE),
           width = 0.4, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0,0), limits = c(0, 1)) +
  scale_fill_manual(
    name = "Cell Type",
    values = setNames(okabe_ito[seq_len(nlevels(ct_prop2$sctype))],
                      levels(ct_prop2$sctype))
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  labs(x = "TB status", y = "cell type proportions (%)", title = NULL) +
  
  # Global font control
  theme_classic(base_size = 10, base_family = "Arial") +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 7),
    legend.key.height = grid::unit(0.35, "cm"),
    legend.key.width  = grid::unit(0.35, "cm"),
    legend.spacing.y  = grid::unit(0.15, "cm"),
    legend.box.margin = margin(0, 0, 0, -10),
    plot.margin = margin(6, 6, 6, 6)
  )
p
# grouped bar side-b side comparison
ggplot(ct_prop2, aes(x = sctype, y = prop, fill = TBstatus)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Case" = "#e6194B", "Control" = "#3cb44b")) +
  labs(x = NULL, y = "Cell type proportion (%)", fill = "TB status") +
  coord_flip() +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(axis.text.y = element_text(face = "bold"))

ggplot(ct_donor, aes(x = sctype, y = prop, fill = TBstatus)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Case" = "#e6194B", "Control" = "#3cb44b")) +
  labs(x = NULL, y = "Cell type proportion (%)", fill = "TB status") +
  coord_flip() +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(axis.text.y = element_text(face = "bold"))


# Now compare distributions, not pooled totals:
ggplot(ct_donor %>% filter(sctype %in% c("Natural killer  cells","CD8+ NKT-like cells")),
       aes(x=TBstatus, y=prop)) +
  geom_boxplot() + geom_jitter(width=0.1, alpha=0.4) +
  scale_y_continuous(labels=scales::percent)




 

DimPlot(seur_obj, group.by = "sctype_classification", label = TRUE, repel = TRUE) +
  NoLegend()
FeaturePlot(seur_obj,
            features = c("NKG7","GNLY","PRF1","KLRD1","FCGR3A","NCR1",
                         "TRAC","CD3D","TRBC1"),
            ncol = 3, order = TRUE)

VlnPlot(seur_obj,
        features = c("TRAC","CD3D","TRBC1","NKG7","GNLY","KLRD1","FCGR3A"),
        group.by = "sctype_classification",
        pt.size = 0, ncol = 2)
DefaultAssay(seur_obj) <- "RNA"

DotPlot(
  seur_obj,
  features = c("TRAC","CD3D","CD247",      # T identity
               "NKG7","GNLY","PRF1",       # cytotoxic program
               "KLRD1","FCGR3A","NCR1"),   # NK identity
  group.by = "sctype_classification"
) + RotatedAxis()



# What assays do you have?
Assays(seur_obj)

# Use RNA feature names
DefaultAssay(seur_obj) <- "RNA"
rna_feats <- rownames(seur_obj[["RNA"]])

# If you have ADT, use those feature names too
adt_feats <- if ("ADT" %in% Assays(seur_obj)) rownames(seur_obj[["ADT"]]) else character(0)

# ScType-style markers (protein-ish + mixed)
markers_raw <- c("CD56","CD16","CD94","NKp46","CD314","CD161","CD122",
                 "NKG7","GNLY","GZMB","GZMA","PRF1","TRAC","CD3D","CD3E","CD3G","CD247")

check =data.frame(
  marker = markers_raw,
  in_RNA = markers_raw %in% rna_feats,
  in_ADT = markers_raw %in% adt_feats
)
VlnPlot(seur_obj,
       features = c("TRAC"),
       group.by = "sctype_classification",
       pt.size = 0, ncol = 1)
#########


