library(dplyr)
library(ggplot2)
library(ggalluvial)

# 1) Make sure these are factors/characters
meta2 <- meta2 %>%
  mutate(across(c(sctype_cd8_nk_tcr_neg, sctype_cd8_nk, sctype_default, sctype_tcr_neg ),
                ~ as.factor(.x)))

# 2) Collapse 500k rows -> unique flows with a count n
meta2_counts <- meta2 %>%
  count(sctype_cd8_nk_tcr_neg, sctype_cd8_nk, sctype_default, sctype_tcr_neg,
        name = "n")

# 3) Plot using y = n (weights)
plt <- ggplot(meta2_counts,
            aes(axis1 = sctype_cd8_nk_tcr_neg,
                axis2 = sctype_cd8_nk,
                axis3 = sctype_default,
                axis4 = sctype_tcr_neg,
                y = n)) +
  geom_alluvium(aes(fill = sctype_default), alpha = 0.8) +
  geom_stratum(color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("sctype_cd8_nk_tcr_neg", "sctype_cd8_nk", "sctype_default", "sctype_tcr_neg"),
    labels = c("CD8/NK (TCR-)", "CD8/NK", "Default", "TCR-"),expand = c(.02, .02)) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Cell Types")

plt <- ggplot(meta2_counts,
            aes(axis1 = sctype_cd8_nk_tcr_neg,
                axis2 = sctype_cd8_nk,
                axis3 = sctype_default,
                axis4 = sctype_default)) +
  geom_alluvium(aes(fill = sctype_default), alpha = 0.8) +
  geom_stratum(color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("sctype_cd8_nk_tcr_neg", "sctype_cd8_nk", "sctype_default", "sctype_tcr_neg"),
    labels = c("CD8/NK (TCR-)", "CD8/NK", "Default", "TCR-"),expand = c(.02, .02)) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(fill = "Cell Types")

ggsave("/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/results/alluvial_wide_vflat.png", plot = plt, width = 20, height = 6.5, units = "in", dpi = 300)
