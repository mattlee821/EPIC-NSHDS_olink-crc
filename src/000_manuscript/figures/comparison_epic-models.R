rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(ggplot2)

# scatter plot comparing results from main analysis with results from datasets of different processing steps ====
## data
data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-FALSE_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt")
data <- data %>%
  dplyr::mutate(
    outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    followup = factor(followup, levels = c(0, 2, 5))
  ) %>%
  dplyr::arrange(exposure, outcome, sex, followup, model) %>%
  dplyr::group_by(exposure, outcome, sex, followup, model) %>%
  dplyr::mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  dplyr::ungroup() 

## data_plot
data_plot <- data %>%
  dplyr::filter(
    followup == 0,
    model %in% c("model_1", "model_2")
  ) %>%
  dplyr::select(exposure, outcome, sex, model, coef, FDR_BH)
## df combined
df1 <- data_plot %>%
  dplyr::filter(model == "model_1") %>%
  dplyr::select(-model)
df2 <- data_plot %>%
  dplyr::filter(model == "model_2") %>%
  dplyr::select(-model)
data_plot <- df1 %>%
  dplyr::rename(coef_data1 = coef, 
         FDR_BH_data1 = FDR_BH) %>%
  dplyr::left_join(
    df2 %>% rename(coef_data2 = coef, 
                   FDR_BH_data2 = FDR_BH),
    by = c("exposure", "outcome", "sex")
  ) %>%
  dplyr::mutate(
    association = case_when(
      FDR_BH_data1 < 0.05 & FDR_BH_data2 < 0.05 ~ "both",   # Both significant
      FDR_BH_data1 < 0.05 ~ "model 1",                     # Only FDR_BH_data1 significant
      FDR_BH_data2 < 0.05 ~ "model 2",                       # Only FDR_BH_data2 significant
      TRUE ~ "none"                                      # All other points
    )
  )

## plot
p1 <- ggplot2::ggplot(data_plot, ggplot2::aes(x = coef_data2, y = coef_data1, colour = association)) +
  ggplot2::geom_point(ggplot2::aes(alpha = ifelse(association == "none", 0.5, 1))) +
  ggplot2::scale_color_manual(
    values = c(
      "none" = "black",    
      "model 1" = "blue",  
      "model 2" = "orange",
      "both" = "green"   
    )) +
  ggplot2::geom_vline(xintercept = 0, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = 0, color = "grey", linetype = "solid") + 
  ggplot2::facet_grid(sex~outcome, 
             scales = "fixed") +
  ggplot2::labs(
    x = "model 1",
    y = "model 2"
  ) +
  cowplot::theme_cowplot() +
  ggplot2::guides(alpha = "none")

grDevices::tiff("manuscript/figures/comparison_models/epic-immuneonc.tiff", 
     width = 800, height = 800, units = "px", res = 100)
print(p1)
grDevices::dev.off()
