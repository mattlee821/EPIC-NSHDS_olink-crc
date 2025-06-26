rm(list=ls())
set.seed(821)

# scatter plot comparing results from main analysis with results from datasets of different processing steps

# Load packages
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)

# List all files
LIST_FILES <- list.files("analysis/002_EPIC-analysis/001_analysis", pattern = ".txt", all.files = T, full.names = T, recursive = T)

# Define the 4 main files
LIST_FILES_main <- c(
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_data-processed_LOD-filter/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LCMD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed_LOD-filter/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LCMD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_intersect_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_intersect_data-processed_LOD-filter/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_intersect_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt",
  "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_intersect_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt"
)

# Build file_info table
file_info <- tibble(file_path = LIST_FILES) %>%
  mutate(
    dataset_type = case_when(
      str_detect(file_path, "EPIC_olink-immuneonc_intersect") ~ "EPIC_olink-immuneonc_intersect",
      str_detect(file_path, "EPIC_olink-explore_intersect") ~ "EPIC_olink-explore_intersect",
      str_detect(file_path, "EPIC_olink-explore") ~ "EPIC_olink-explore",
      str_detect(file_path, "EPIC_olink-immuneonc") ~ "EPIC_olink-immuneonc",
      TRUE ~ NA_character_
    ),
    processing_type = case_when(
      str_detect(file_path, "data-processed_LOD-filter") ~ "data-processed_LOD-filter",
      str_detect(file_path, "data-processed/") ~ "data-processed",
      TRUE ~ NA_character_
    ),
    file_category = if_else(file_path %in% LIST_FILES_main, "main", "other")
  ) %>%
  filter(!is.na(dataset_type), !is.na(processing_type))

# Initialize result table
table <- data.frame(outcome = character(), sex = character(), cor_all = numeric(), cor_fdr = numeric())

# Loop through main files
for (main_file in LIST_FILES_main) {
  # Inside the for-loop for each main file
  main_info <- file_info %>% filter(file_path == main_file)
  main_dataset_type <- main_info$dataset_type
  main_processing_type <- main_info$processing_type
  
  # Correct filtering of other files
  other_files <- file_info %>%
    filter(
      dataset_type == main_dataset_type,
      processing_type == main_processing_type,
      file_category == "other"
    ) %>%
    pull(file_path)
  
  # Label for main file
  label <- basename(main_file) %>%
    gsub(".txt$", "", .) %>%
    gsub("data-features_", "", .)
  
  # Read in main data
  data_main <- fread(main_file) %>%
    mutate(
      outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
      sex = factor(sex, levels = c("combined", "female", "male")),
      followup = factor(followup, levels = c(0, 2, 5))
    ) %>%
    arrange(exposure, outcome, sex, followup, model) %>%
    group_by(exposure, outcome, sex, followup, model) %>%
    mutate(FDR_BH = p.adjust(pval, method = "BH")) %>%
    ungroup() %>%
    filter(followup == 0, model == "model_2")
  
  # Loop through matched other files
  for (other_file in other_files) {
    cat("Comparing main with:", basename(other_file), "\n")
    
    label_other <- basename(other_file) %>%
      gsub(".txt$", "", .) %>%
      gsub("data-features_", "", .)
    
    data_other <- fread(other_file) %>%
      mutate(
        outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
        sex = factor(sex, levels = c("combined", "female", "male")),
        followup = factor(followup, levels = c(0, 2, 5))
      ) %>%
      arrange(exposure, outcome, sex, followup, model) %>%
      group_by(exposure, outcome, sex, followup, model) %>%
      mutate(FDR_BH = p.adjust(pval, method = "BH")) %>%
      ungroup() %>%
      filter(followup == 0, model == "model_2")
    
    # Combine datasets
    data_plot <- data_main %>%
      select(exposure, outcome, sex, coef_data1 = coef, FDR_BH_data1 = FDR_BH) %>%
      inner_join(data_other %>%
                   select(exposure, outcome, sex, coef_data2 = coef, FDR_BH_data2 = FDR_BH),
                 by = c("exposure", "outcome", "sex")) %>%
      mutate(
        association = case_when(
          FDR_BH_data1 < 0.05 & FDR_BH_data2 < 0.05 ~ "both",
          FDR_BH_data1 < 0.05 ~ "main",
          FDR_BH_data2 < 0.05 ~ "sensitivity",
          TRUE ~ "none"
        )
      )
    
    # Correlation summary
    table_cor <- data_plot %>%
      group_by(outcome, sex) %>%
      summarise(
        cor_all = cor(coef_data1, coef_data2, use = "complete.obs", method = "pearson"),
        cor_fdr = cor(coef_data1[FDR_BH_data1 < 0.05], coef_data2[FDR_BH_data1 < 0.05], use = "complete.obs", method = "pearson"),
        n_proteins = n_distinct(exposure),
        .groups = "drop"
      ) %>%
      mutate(dataset = label, 
             dataset_type = main_dataset_type,
             processing_type = main_processing_type) %>%
      select(dataset, processing_type, outcome, sex, cor_all, cor_fdr, n_proteins)
    
    main_info <- file_info %>% filter(file_path == main_file)
    main_dataset_type <- main_info$dataset_type
    main_processing_type <- main_info$processing_type
    
    table <- bind_rows(table, table_cor)
    
    # Plotting
    out_dir <- file.path("manuscript", "figures", "comparison_datasets", main_dataset_type, main_processing_type)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_dir, paste0(label_other, ".pdf"))
    
    pdf(out_file, width = 14, height = 14)  # open a PDF file
    for (sex_val in unique(data_plot$sex)) {
      for (outcome_val in unique(data_plot$outcome)) {
        
        sub_data <- data_plot %>%
          dplyr::filter(sex == sex_val, outcome == outcome_val)
        
        cor_row <- table_cor %>%
          dplyr::filter(sex == sex_val, outcome == outcome_val)
        
        title_text <- paste0(
          "Sex: ", sex_val,
          ", Outcome: ", outcome_val,
          "\nN proteins: ", cor_row$n_proteins,
          ", cor_all: ", round(cor_row$cor_all, 3),
          ", cor_fdr: ", round(cor_row$cor_fdr, 3)
        )
        
        p <- ggplot2::ggplot(sub_data, ggplot2::aes(x = coef_data1, y = coef_data2, colour = association)) +
          ggplot2::geom_point(ggplot2::aes(alpha = ifelse(association == "none", 0.5, 1))) +
          ggplot2::scale_color_manual(values = c("both" = "green", "main" = "orange", "sensitivity" = "blue", "none" = "black")) +
          ggplot2::geom_vline(xintercept = 0, color = "grey") +
          ggplot2::geom_hline(yintercept = 0, color = "grey") +
          ggplot2::labs(title = title_text, x = label_other, y = label) +
          cowplot::theme_cowplot() +
          ggplot2::guides(alpha = "none") 
        print(p)
      }
    }
    dev.off()
  }
}

write.table(table, "manuscript/tables/comparison_epic-datasets_correlations.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
