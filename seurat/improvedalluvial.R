library(dplyr)
library(ggplot2)
library(ggalluvial)

meta =seur_obj@meta.data
# 1) Make sure these are factors/characters
meta2 <- meta %>%
  mutate(across(c(sctype_cd8_nk_tcr_neg, sctype_cd8_nk, sctype_default, sctype_tcr_neg , predicted.celltype.l2),
                ~ as.factor(.x)))

# 2) Collapse 500k rows -> unique flows with a count n
meta2_counts <- meta2 %>%
  count(sctype_cd8_nk_tcr_neg, sctype_cd8_nk, sctype_default, sctype_tcr_neg,sctype_oscar, predicted.celltype.l2,
        name = "n")  

# meta2_counts <- meta2 %>%
#   count(sctype_cd8_nk_tcr_neg, sctype_cd8_nk, sctype_default,
#         sctype_tcr_neg, sctype_oscar, name = "n") %>%
#   filter(!sctype_default %in% c(
#     "Non-classical monocytes","Plasma B cells","Plasmacytoid Dendritic cells",
#     "Memory CD4+ T cells","Naive CD4+ T cells","Myeloid Dendritic cells", "Classical Monocytes", "Naive B cells","Effector CD8+ T cells",
#     "Megakaryocyte", "Naive CD8+ T cells", "Effector CD8+ T cells"
#   ))

# alluvial_cell_plt <- ggplot(meta2_counts,
#                             aes(axis1 = sctype_default,
#                                 axis2 = predicted.celltype.l2)) +
#   geom_alluvium(aes(fill = sctype_default), alpha = 0.8) +
#   geom_stratum(color = "grey30") +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
#   scale_x_discrete(
#     limits = c( "sctype_default", "predicted.celltype.l2"),
#     labels = c("Default", "Azimuth"),expand = c(.02, .02)) +
#   theme_bw(base_size = 12) +
#   theme(axis.title = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()) +
#   labs(fill = "Cell Types - default")

# 3) Plot using y = n (weights)
 

alluvial_cell_plt <- ggplot(meta2_counts,
            aes(axis1 = sctype_cd8_nk_tcr_neg,
                axis2 = sctype_cd8_nk,
                axis3 = sctype_default,
                axis4 = sctype_tcr_neg,
                axis5 = sctype_oscar)) +
  geom_alluvium(aes(fill = sctype_default), alpha = 0.8) +
  geom_stratum(color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("sctype_cd8_nk_tcr_neg", "sctype_cd8_nk", "sctype_default", "sctype_tcr_neg","sctype_oscar"),
    labels = c("CD8/NK (TCR-)", "CD8/NK", "Default", "TCR-", "Oscar"),expand = c(.02, .02)) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Cell Types - default")


 


ggsave("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/alluvial_pbmcs.png", plot = alluvial_cell_plt, width = 20, height = 6.5, units = "in", dpi = 300)
#ggsave("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/alluvial_pbmcsunflt.png", plot = alluvial_cell_plt, width = 20, height = 6.5, units = "in", dpi = 300)



make_alluvial_one <- function(df_one, label, theme_obj = NULL) {
  p <- ggplot(df_one,
              aes(axis1 = sctype_default,
                  axis2 = predicted.celltype.l2,
                  y = n)) +
    geom_alluvium(aes(fill = sctype_default), alpha = 0.8) +
    geom_stratum(color = "grey30") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(
      limits = c("sctype_default", "predicted.celltype.l2"),
      labels = c("ScType", "Azimuth"),
      expand = c(.02, .02)
    ) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = paste0("ScType: ", label),
         fill  = "ScType label")
  
  if (!is.null(theme_obj)) p <- p + theme_obj
  p
}

sctype_levels <- sort(unique(as.character(meta2_counts$sctype_default)))

plots_by_sctype <- lapply(sctype_levels, function(lbl) {
  df_one <- meta2_counts %>% filter(sctype_default == lbl)
  make_alluvial_one(df_one, label = lbl, theme_obj = NULL)  # or NULL
})
names(plots_by_sctype) <- sctype_levels

saveRDS(plots_by_sctype,"/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/AZIMUTHPLOTS.rds")



