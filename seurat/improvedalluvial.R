library(dplyr)
library(ggplot2)
library(ggalluvial)

# 1) Make sure these are factors/characters
meta2 <- meta %>%
  mutate(across(c(sctype_cd8_nk_tcr_neg, sctype_cd8_nk, sctype_default, sctype_tcr_neg ),
                ~ as.factor(.x)))

# 2) Collapse 500k rows -> unique flows with a count n
meta2_counts <- meta2 %>%
  count(sctype_cd8_nk_tcr_neg, sctype_cd8_nk, sctype_default, sctype_tcr_neg,sctype_oscar,
        name = "n")  

meta2_counts <- meta2 %>%
  count(sctype_cd8_nk_tcr_neg, sctype_cd8_nk, sctype_default,
        sctype_tcr_neg, sctype_oscar, name = "n") %>%
  filter(!sctype_default %in% c(
    "Non-classical monocytes","Plasma B cells","Plasmacytoid Dendritic cells",
    "Memory CD4+ T cells","Naive CD4+ T cells","Myeloid Dendritic cells", "Classical Monocytes", "Naive B cells", "γδ-T cells","Unknown" ,"Effector CD8+ T cells",
    "Megakaryocyte", "Naive CD8+ T cells", "Effector CD8+ T cells"
  ))

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
