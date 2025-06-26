rm(list=ls())
set.seed(821)

# environment ====
source("src/000_source.R")
rm(list=ls())
source("src/000_functions.R")
cat("\n # START")

# data ====
data_immuneonc <- data.table::fread("data/processed/EPIC_olink-immuneonc.txt") %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  tibble::as_tibble() 
ID_immuneonc <- data_immuneonc %>%
  dplyr::select(SampleID) %>%
  unique()
feature_immuneonc <- data_immuneonc %>%
  dplyr::select(UniProt) %>%
  unique()

data_explore <- data.table::fread("data/processed/EPIC_olink-explore.txt") %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  tibble::as_tibble() 
ID_explore <- data_explore %>%
  dplyr::select(SampleID) %>%
  unique()
feature_explore <- data_explore %>%
  dplyr::select(UniProt) %>%
  unique()

ID_sample <- intersect(ID_explore, ID_immuneonc) %>%
  dplyr::pull(SampleID)
ID_feature <- intersect(feature_immuneonc, feature_explore) %>%
  dplyr::pull(UniProt)

# filter for intersect ====
data_immuneonc <- data_immuneonc %>%
  dplyr::filter(SampleID %in% ID_sample) %>%
  dplyr::filter(UniProt %in% ID_feature)

data_explore <- data_explore %>%
  dplyr::filter(SampleID %in% ID_sample) %>%
  dplyr::filter(UniProt %in% ID_feature)

# data_immuneonc ====
data <- data_immuneonc
data_features <- data %>%
  dplyr::select(SampleID, OlinkID, NPX) %>%
  unique() %>%
  tidyr::pivot_wider(names_from = OlinkID, values_from = NPX)

data_meta_samples <- data %>%
  dplyr::select(-c(2:14)) %>%
  unique

data_meta_features <- data.table::fread("analysis/001_data-processing/EPIC_olink-immuneonc_metadata-within-batches.txt") %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  tibble::as_tibble() %>%
  dplyr::select(SampleID, feature, LOD) %>%
  unique()

# identify features with > X% values < LOD
features_below_LOD <- identify_features_below_lod(
  data_features = data_features, 
  data_meta_features = data_meta_features,
  feature_col = "feature", 
  LOD_col = "LOD",
  percent = 10)

df_below_lod_percentages_features <- calculate_below_lod_percentages(
  data_features = data_features, 
  data_meta_features = data_meta_features,
  feature_col = "feature", 
  LOD_col = "LOD")

p1 <- plot_features_exceeding_threshold(
  df_below_lod_percentages = df_below_lod_percentages_features)

### samples
samples_below_LOD <- identify_samples_below_lod(
  data_features = data_features, 
  data_meta_features = data_meta_features,
  feature_col = "feature", 
  sample_col = "SampleID",
  LOD_col = "LOD", 
  percent = 10)

df_below_lod_percentages_samples <- calculate_below_lod_percentages_samples(
  data_features = data_features,
  data_meta_features = data_meta_features,
  feature_col = "feature",
  sample_col = "SampleID",
  LOD_col = "LOD")

p2 <- plot_samples_exceeding_threshold(
  df_below_lod_percentages = df_below_lod_percentages_samples, 
  sample_col = "SampleID")

### combined plot 
grDevices::tiff(filename = "analysis/001_data-processing/EPIC_olink-immuneonc_explore-intersect_barchart-LOD.tiff", 
                width = 800, height = 600, units = "px")
cowplot::plot_grid(p1,p2, nrow = 1)
grDevices::dev.off()

# create df with LOD filter
data_features_LOD_filter <- data_features %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "feature",
    values_to = "value"
  ) %>%
  left_join(data_meta_features, by = c("SampleID", "feature")) %>%
  mutate(value = ifelse(value < LOD, NA, value)) %>%
  select(SampleID, feature, value) %>%
  pivot_wider(names_from = feature, values_from = value) %>%
  unique()

# process data
data_list <- list(
  "data_features" = list(data = data_features, dir = "data-processed"),
  "data_features_LOD_filter" = list(data = data_features_LOD_filter, dir = "data-processed_LOD-filter")
)

for (name in names(data_list)) {
  df <- data_list[[name]]
  
  for (transform_flag in c(TRUE, FALSE)) {
    for (outlier_flag in c(TRUE, FALSE)) {
      for (imputation_flag in c(TRUE, FALSE)) {
        for (imputation_method_flag in c("LOD", "LCMD")) {
          
        DIR <- paste0("analysis/001_data-processing/EPIC_olink-immuneonc_intersect_", df$dir, "/")
        
        cat("\n ############# \n", sep = "")
        
        OmicsProcessing::process_data(
          data = df$data, 
          data_meta_features = data_meta_features, 
          data_meta_samples = data_meta_samples, 
          col_samples = "SampleID", 
          col_features = "feature", 
          
          save = TRUE, 
          path_out = DIR, 
          path_outliers = DIR,
          
          # Constants
          exclusion_extreme_feature = TRUE, 
          missing_pct_feature = 0.8, 
          exclusion_extreme_sample = TRUE, 
          missing_pct_sample = 0.8, 
          
          imputation = imputation_flag, 
          imputation_method = imputation_method_flag, 
          col_LOD = "LOD", 
          
          transformation = transform_flag, 
          transformation_method = "InvRank", 
          
          outlier = outlier_flag, 
          
          plate_correction = FALSE, 
          cols_listRandom = NULL, 
          cols_listFixedToKeep = NULL, 
          cols_listFixedToRemove = NULL, 
          col_HeteroSked = NULL, 
          
          centre_scale = TRUE, 
          
          case_control = TRUE, 
          col_case_control = "Match_Caseset"
        )
        cat("\n")
      }
    }
  }
  }
}

# data_explore ====
data <- data_explore
data_features <- data %>%
  dplyr::select(SampleID, OlinkID, NPX) %>%
  unique() %>%
  tidyr::pivot_wider(names_from = OlinkID, values_from = NPX)

data_meta_samples <- data %>%
  dplyr::select(-c(2:16)) %>%
  unique

data_meta_features <- data.table::fread("analysis/001_data-processing/EPIC_olink-explore_metadata-within-batches.txt") %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  tibble::as_tibble() %>%
  dplyr::select(SampleID, feature, LOD) %>%
  unique()

# identify features with > X% values < LOD
features_below_LOD <- identify_features_below_lod(
  data_features = data_features, 
  data_meta_features = data_meta_features,
  feature_col = "feature", 
  LOD_col = "LOD",
  percent = 10)

df_below_lod_percentages_features <- calculate_below_lod_percentages(
  data_features = data_features, 
  data_meta_features = data_meta_features,
  feature_col = "feature", 
  LOD_col = "LOD")

p1 <- plot_features_exceeding_threshold(
  df_below_lod_percentages = df_below_lod_percentages_features)

### samples
samples_below_LOD <- identify_samples_below_lod(
  data_features = data_features, 
  data_meta_features = data_meta_features,
  feature_col = "feature", 
  sample_col = "SampleID",
  LOD_col = "LOD", 
  percent = 10)

df_below_lod_percentages_samples <- calculate_below_lod_percentages_samples(
  data_features = data_features,
  data_meta_features = data_meta_features,
  feature_col = "feature",
  sample_col = "SampleID",
  LOD_col = "LOD")

p2 <- plot_samples_exceeding_threshold(
  df_below_lod_percentages = df_below_lod_percentages_samples, 
  sample_col = "SampleID")

### combined plot 
grDevices::tiff(filename = "analysis/001_data-processing/EPIC_olink-explore_immuneonc-intersect_barchart-LOD.tiff", 
                width = 800, height = 600, units = "px")
cowplot::plot_grid(p1,p2, nrow = 1)
grDevices::dev.off()

# create df with LOD filter
data_features_LOD_filter <- data_features %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "feature",
    values_to = "value"
  ) %>%
  left_join(data_meta_features, by = c("SampleID", "feature")) %>%
  mutate(value = ifelse(value < LOD, NA, value)) %>%
  select(SampleID, feature, value) %>%
  pivot_wider(names_from = feature, values_from = value) %>%
  unique()

# process data
data_list <- list(
  "data_features" = list(data = data_features, dir = "data-processed"),
  "data_features_LOD_filter" = list(data = data_features_LOD_filter, dir = "data-processed_LOD-filter")
)

for (name in names(data_list)) {
  df <- data_list[[name]]
  
  for (transform_flag in c(TRUE, FALSE)) {
    for (outlier_flag in c(TRUE, FALSE)) {
      for (imputation_flag in c(TRUE, FALSE)) {
        for (imputation_method_flag in c("LOD", "LCMD")) {
          
          DIR <- paste0("analysis/001_data-processing/EPIC_olink-explore_intersect_", df$dir, "/")
          
          cat("\n ############# \n", sep = "")
          
          OmicsProcessing::process_data(
            data = df$data, 
            data_meta_features = data_meta_features, 
            data_meta_samples = data_meta_samples, 
            col_samples = "SampleID", 
            col_features = "feature", 
            
            save = TRUE, 
            path_out = DIR, 
            path_outliers = DIR,
            
            # Constants
            exclusion_extreme_feature = TRUE, 
            missing_pct_feature = 0.8, 
            exclusion_extreme_sample = TRUE, 
            missing_pct_sample = 0.8, 
            
            imputation = imputation_flag, 
            imputation_method = imputation_method_flag, 
            col_LOD = "LOD", 
            
            transformation = transform_flag, 
            transformation_method = "InvRank", 
            
            outlier = outlier_flag, 
            
            plate_correction = FALSE, 
            cols_listRandom = NULL, 
            cols_listFixedToKeep = NULL, 
            cols_listFixedToRemove = NULL, 
            col_HeteroSked = NULL, 
            
            centre_scale = TRUE, 
            
            case_control = TRUE, 
            col_case_control = "Match_Caseset"
          )
          cat("\n")
        }
      }
    }
  }
}

# ====
cat("\n # END")
