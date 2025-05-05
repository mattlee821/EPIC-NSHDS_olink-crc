rm(list=ls())
set.seed(821)

# environment ====
source("src/000_functions.R")

# data ====
data_meta_features <- data.table::fread("analysis/001_data-processing/EPIC_olink-immuneonc_feature-metadata.txt")
data_meta_samples <- data.table::fread("analysis/001_data-processing/EPIC_olink-immuneonc_sample-metadata.txt")
data_features <- data.table::fread("data/processed/EPIC_olink-immuneonc_feature-data.txt")

# processing:
# 1. process as is
# 2. process after replacing all values below LOD with the LOD
# 3. process after replacing all values below LOD with NA

# imputation test ====
## to identify which imputation method performs best for our data
## based on this report RF looks best
ImputationReport::imputation_test(
  data = data_features, 
  knit_root_dir = "/data/IARC_Biostat/work/EPIC-NSHDS_olink-crc/", 
  output_dir = "analysis/001_data-processing/",
  subtitle = "EPIC"
)

# identify features with > X% values < LOD ====
## the below shows that there are 10
### features
features_below_LOD <- identify_features_below_lod(
  data_features = data_features, 
  data_meta_features = data_meta_features,
  feature_col = "OlinkID", 
  LOD_col = "LOD",
  percent = 10)

df_below_lod_percentages_features <- calculate_below_lod_percentages(
  data_features = data_features, 
  data_meta_features = data_meta_features, 
  feature_col = "OlinkID", 
  LOD_col = "LOD")

p1 <- plot_features_exceeding_threshold(
  df_below_lod_percentages = df_below_lod_percentages_features)

### samples
samples_below_LOD <- identify_samples_below_lod(
  data_features = data_features, 
  data_meta_features = data_meta_features, 
  feature_col = "OlinkID", 
  sample_col = "Idepic_Bio",
  LOD_col = "LOD", 
  percent = 10)

df_below_lod_percentages_samples <- calculate_below_lod_percentages_samples(
  data_features = data_features,
  data_meta_features = data_meta_features,
  feature_col = "OlinkID",
  LOD_col = "LOD",
  sample_col = "Idepic_Bio")

p2 <- plot_samples_exceeding_threshold(
  df_below_lod_percentages = df_below_lod_percentages_samples, 
  sample_col = "Idepic_Bio")

### combined plot 
grDevices::tiff(filename = "analysis/001_data-processing/EPIC_olink-immuneonc_barchart-LOD.tiff", 
                width = 800, height = 600, units = "px")
cowplot::plot_grid(p1,p2, nrow = 1)
grDevices::dev.off()

# create df with LOD filter ====
data_features_LOD_filter <- data_features %>%
  dplyr::select(!dplyr::all_of(
    df_below_lod_percentages_features %>%
      dplyr::filter(Percentage_Below_LOD > 40) %>% # remove features with >X% of values < LOD (using 40% as its nested case-control so 50% could correspond to case-status)
      dplyr::pull(Feature)
  ))

# process data ====
## needs to be done for data_features and data_features_LOD_filter
temp <- OmicsProcessing::process_data(
  data = data_features, 
  data_meta_features = data_meta_features, 
  data_meta_samples = data_meta_samples, 
  col_samples = "Idepic_Bio", 
  col_features = "OlinkID", 
  
  save = TRUE, 
  path_out = "analysis/001_data-processing/EPIC_olink-immuneonc_data-processed/", 
  path_outliers = "analysis/001_data-processing/EPIC_olink-immuneonc_data-processed/", 
  
  ## dont change
  exclusion_extreme_feature = TRUE, 
  missing_pct_feature = 0.8, 
  exclusion_extreme_sample = TRUE, 
  missing_pct_sample = 0.8, 
  
  ## dont change
  imputation = FALSE, 
  imputation_method = "LOD", 
  col_LOD = "LOD", 
  
  ## do no transformation and InvRank transformation
  transformation = TRUE, 
  transformation_method = "InvRank", 
  
  ## dont change
  outlier = TRUE, 
  
  ## dont change
  ## based on analysis/001_data-processing/EPIC_olink-immuneonc_scatter-batch.tiff there appears to be no batch effects so we wont do plate_correction
  plate_correction = FALSE, 
  cols_listRandom = NULL, 
  cols_listFixedToKeep = NULL, 
  cols_listFixedToRemove = NULL, 
  col_HeteroSked = NULL, 
  
  ## dont change
  centre_scale = TRUE, 
  
  ## dont change
  case_control = TRUE, 
  col_case_control = "Match_Caseset")
