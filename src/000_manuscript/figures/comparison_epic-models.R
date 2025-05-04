rm(list=ls())
set.seed(821)

# environment ====
library(functions)
library(ggplot2)
library(wesanderson)
library(dplyr)
library(cowplot)
library(data.table)
library(ggrepel)
library(tidyr)
library(data.table)
library(dplyr)
library(tidyr)
library(pathfindR)
library(cowplot)
library(grid)
library(UpSetR)
library(ggplot2)

# source ====
palette <- palette()

# scatter plot comparing results from main analysis with results from datasets of different processing steps ====
## data
data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-FALSE_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt")
data <- data %>%
  mutate(
    outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    followup = factor(followup, levels = c(0, 2, 5))
  ) %>%
  arrange(exposure, outcome, sex, followup, model) %>%
  group_by(exposure, outcome, sex, followup, model) %>%
  mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  ungroup() 

## data_plot
data_plot <- data %>%
  filter(
    followup == 0,
    model %in% c("model_1", "model_2")
  ) %>%
  select(exposure, outcome, sex, model, coef, FDR_BH)
## df combined
df1 <- data_plot %>%
  filter(model == "model_1") %>%
  select(-model)
df2 <- data_plot %>%
  filter(model == "model_2") %>%
  select(-model)
data_plot <- df1 %>%
  rename(coef_data1 = coef, 
         FDR_BH_data1 = FDR_BH) %>%
  left_join(
    df2 %>% rename(coef_data2 = coef, 
                   FDR_BH_data2 = FDR_BH),
    by = c("exposure", "outcome", "sex")
  ) %>%
  mutate(
    association = case_when(
      FDR_BH_data1 < 0.05 & FDR_BH_data2 < 0.05 ~ "both",   # Both significant
      FDR_BH_data1 < 0.05 ~ "model 1",                     # Only FDR_BH_data1 significant
      FDR_BH_data2 < 0.05 ~ "model 2",                       # Only FDR_BH_data2 significant
      TRUE ~ "none"                                      # All other points
    )
  )

## plot
p1 <- ggplot(data_plot, aes(x = coef_data2, y = coef_data1, colour = association)) +
  geom_point(aes(alpha = ifelse(association == "none", 0.5, 1))) +
  scale_color_manual(
    values = c(
      "none" = "black",    
      "model 1" = "blue",  
      "model 2" = "orange",
      "both" = "green"   
    )) +
  geom_vline(xintercept = 0, color = palette[[1]][[1]], linetype = "solid") +
  geom_hline(yintercept = 0, color = palette[[1]][[1]], linetype = "solid") + 
  facet_grid(sex~outcome, 
             scales = "fixed") +
  labs(
    x = "model 1",
    y = "model 2"
  ) +
  cowplot::theme_cowplot() +
  guides(alpha = "none")

tiff("manuscript/figures/comparison_models/epic-immuneonc.tiff", 
     width = 800, height = 800, units = "px", res = 100)
print(p1)
dev.off()

