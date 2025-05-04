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
list_files <- list.files(path = "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed", pattern = "data-features", all.files = T, recursive = T, full.names = T)
main_file <- "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-FALSE_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt"
data <- data.table::fread(main_file)
data <- data %>%
  mutate(
    outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    followup = factor(followup, levels = c(0, 2, 5))
  ) %>%
  arrange(exposure, outcome, sex, followup, model) %>%
  group_by(exposure, outcome, sex, followup, model) %>%
  mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  ungroup() %>%
  filter(followup == 0, model == "model_2")
## table/plot loop
table <- data.frame(data = character(), outcome = character(), sex = character(), cor_all = numeric(), cor_fdr = numeric())

# Process each file sequentially
for (file in list_files) {
  # Skip the main file
  if (file == main_file) next
  label <- basename(file)
  label <- gsub(".txt", "", label)
  label <- gsub("data-features_", "", label)
  # Read and format the new data
  df <- data.table::fread(file) 
  df <- df %>%
    mutate(
      outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
      sex = factor(sex, levels = c("combined", "female", "male")),
      followup = factor(followup, levels = c(0, 2, 5))
    ) %>%
    arrange(exposure, outcome, sex, followup, model) %>%
    group_by(exposure, outcome, sex, followup, model) %>%
    mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
    ungroup() %>%
    filter(followup == 0, model == "model_2")
  
  # combined
  data_plot <- data %>%
    select(exposure, outcome, sex, coef_data1 = coef, FDR_BH_data1 = FDR_BH) %>%
    inner_join(
      df %>% select(exposure, outcome, sex, coef_data2 = coef, FDR_BH_data2 = FDR_BH),
      by = c("exposure", "outcome", "sex")
    ) %>%
    mutate(
      association = case_when(
        FDR_BH_data1 < 0.05 & FDR_BH_data2 < 0.05 ~ "both",   # Both significant
        FDR_BH_data1 < 0.05 ~ "main",                     # Only FDR_BH_data1 significant
        FDR_BH_data2 < 0.05 ~ "sensitivity",                       # Only FDR_BH_data2 significant
        TRUE ~ "none"                                      # All other points
      )
    )
  
  # Compute correlation between coef columns
  table_cor <- data_plot %>%
    group_by(outcome, sex) %>%
    summarise(
      cor_all = cor(coef_data1, coef_data2, use = "complete.obs", method = "pearson"),  # Full data correlation
      cor_fdr = cor(coef_data1[FDR_BH_data1 < 0.05], coef_data2[FDR_BH_data1 < 0.05], use = "complete.obs", method = "pearson"),  # Significant data correlation
      .groups = "drop"
    ) %>%
    mutate(dataset = label) %>%
    select(dataset, outcome, sex, cor_all, cor_fdr)
  table <- bind_rows(table, table_cor)
  cat("cor done \n")
  
  # plot
  p1 <- ggplot(data_plot, aes(x = coef_data2, y = coef_data1, colour = association)) +
    geom_point(aes(alpha = ifelse(association == "none", 0.5, 1))) +  
    scale_color_manual(
      values = c(
        "both" = "green",   
        "main" = "orange",  
        "sensitivity" = "blue",
        "none" = "black"    
      )) +
    geom_vline(xintercept = 0, color = palette[[1]][[1]], linetype = "solid") +
    geom_hline(yintercept = 0, color = palette[[1]][[1]], linetype = "solid") + 
    facet_grid(sex~outcome, 
               scales = "fixed") +
    labs(
      x = label,
      y = "Main analysis data"
    ) +
    cowplot::theme_cowplot() +
    guides(alpha = "none")  
  tiff(paste0("manuscript/figures/comparison_datasets/epic-immuneonc/", label, ".tiff"), 
       width = 1400, height = 1000, units = "px", res = 100)
  print(p1)
  dev.off()
  cat("plot done \n")
  
}
