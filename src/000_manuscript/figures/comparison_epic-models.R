rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(data.table)

# files ====
LIST_FILES <- c(
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_data-processed_LOD-filter/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LCMD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed_LOD-filter/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LCMD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt"
)
file_info <- tibble(
  file_path = LIST_FILES
) %>%
  mutate(
    dataset_type = case_when(
      str_detect(file_path, "EPIC_olink-explore") ~ "EPIC_olink-explore",
      str_detect(file_path, "EPIC_olink-immuneonc") ~ "EPIC_olink-immuneonc",
      TRUE ~ NA_character_
    ),
    processing_type = case_when(
      str_detect(file_path, "data-processed_LOD-filter") ~ "data-processed_LOD-filter",
      str_detect(file_path, "data-processed/") ~ "data-processed",
      TRUE ~ NA_character_
    ))

# function ====
process_file <- function(file_path, dataset_type, processing_type) {
  cat("Processing file:", basename(file_path), "\n")
  
  data <- fread(file_path) %>%
    mutate(
      outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
      sex = factor(sex, levels = c("combined", "female", "male")),
      followup = factor(followup, levels = c(0, 2, 5))
    ) %>%
    arrange(exposure, outcome, sex, followup, model) %>%
    group_by(exposure, outcome, sex, followup, model) %>%
    mutate(FDR_BH = p.adjust(pval, method = "BH")) %>%
    ungroup()

    # Filter model_1 and model_2 at followup == 0
    data_plot <- data %>%
      filter(followup == 0, model %in% c("model_1", "model_2")) %>%
      select(exposure, outcome, sex, model, coef, FDR_BH)
    
    df1 <- data_plot %>% filter(model == "model_1") %>% select(-model)
    df2 <- data_plot %>% filter(model == "model_2") %>% select(-model)
    
    combined <- df1 %>%
      rename(coef_model1 = coef, FDR_BH_model1 = FDR_BH) %>%
      left_join(
        df2 %>% rename(coef_model2 = coef, FDR_BH_model2 = FDR_BH),
        by = c("exposure", "outcome", "sex")
      ) %>%
      mutate(
        association = case_when(
          FDR_BH_model1 < 0.05 & FDR_BH_model2 < 0.05 ~ "both",
          FDR_BH_model1 < 0.05 ~ "model 1",
          FDR_BH_model2 < 0.05 ~ "model 2",
          TRUE ~ "none"
        )
      )
    
    unique_sex <- unique(combined$sex)
    unique_outcome <- unique(combined$outcome)
    
    for (sex_val in unique_sex) {
      for (outcome_val in unique_outcome) {
        subset_data <- combined %>% filter(sex == sex_val, outcome == outcome_val)
        
        cor_all <- cor(subset_data$coef_model1, subset_data$coef_model2, use = "complete.obs")
        cor_fdr <- cor(
          subset_data$coef_model1[subset_data$FDR_BH_model1 < 0.05],
          subset_data$coef_model2[subset_data$FDR_BH_model1 < 0.05],
          use = "complete.obs"
        )
        n_proteins <- n_distinct(subset_data$exposure)
        
        title_text <- paste0(
          "Dataset: ", dataset_type, ", Processing: ", processing_type, "\n",
          "Sex: ", sex_val, ", Outcome: ", outcome_val, "\n",
          "N proteins: ", n_proteins, 
          ", cor_all: ", round(cor_all, 3), 
          ", cor_fdr: ", round(cor_fdr, 3)
        )
        
        p <- ggplot(subset_data, aes(x = coef_model2, y = coef_model1, colour = association)) +
          geom_point(aes(alpha = ifelse(association == "none", 0.5, 1))) +
          scale_color_manual(values = c("both" = "green", "model 1" = "blue", "model 2" = "orange", "none" = "black")) +
          geom_vline(xintercept = 0, color = "grey") +
          geom_hline(yintercept = 0, color = "grey") +
          labs(title = title_text, x = "model 2", y = "model 1") +
          cowplot::theme_cowplot() +
          guides(alpha = "none")
        
        out_dir <- file.path("manuscript", "figures", "comparison_models", dataset_type, processing_type)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        
        out_file <- file.path(out_dir, paste0(sex_val, "_", outcome_val, ".tiff"))
        
        tiff(out_file, width = 1400, height = 1000, units = "px", res = 100)
        print(p)
        dev.off()
        
        cat("    Plot saved for sex =", sex_val, "and outcome =", outcome_val, "\n")
      }
    }
}

# make plots ====
file_info %>% 
  filter(!is.na(dataset_type) & !is.na(processing_type)) %>% 
  rowwise() %>% 
  do({
    process_file(.$file_path, .$dataset_type, .$processing_type)
    data.frame()
  })
