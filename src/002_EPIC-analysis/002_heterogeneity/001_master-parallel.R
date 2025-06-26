rm(list=ls())
set.seed(821)

#environment ====
source("src/000_source.R")
rm(list=ls())
source("src/000_functions.R")

# Files ====
LIST_FILES <- list.files("analysis/001_data-processing/", pattern = "data-features_", all.files = T, full.names = T, recursive = T)

# Analysis function ====
process_file <- function(i) {
  VAR <- stringr::str_match(
    i,
    ".*/(EPIC_olink-[^/]+)_(data-processed(?:_LOD-filter)?)/.*"
  )
  
  # Data ====
  cat("Processing: \n", basename(i), "\n", sep = "")
  load(i)
  
  data_features <- df_processed
  
  data_meta_samples <- data.table::fread(paste0("data/processed/", gsub(pattern = "_intersect", replacement = "", x = VAR[2]), ".txt")) %>%
    dplyr::filter(Sample_Type == "SAMPLE") %>%
    tibble::as_tibble() %>%
    dplyr::filter(SampleID %in% data_features$SampleID) 
  # Conditional column selection
  if (grepl("explore", VAR[2], ignore.case = TRUE)) {
    data_meta_samples <- data_meta_samples %>%
      dplyr::select(-c(2:16)) %>%
      unique()
  } else if (grepl("immuneonc", VAR[2], ignore.case = TRUE)) {
    data_meta_samples <- data_meta_samples %>%
      dplyr::select(-c(2:14)) %>%
      unique()
  } else {
    # Fallback option (optional)
    data_meta_samples <- data_meta_samples %>%
      unique()
  }
  
  data_meta_features <- data.table::fread(paste0("analysis/001_data-processing/", gsub(pattern = "_intersect", replacement = "", x = VAR[2]), "_metadata-within-batches.txt")) %>%
    dplyr::filter(Sample_Type == "SAMPLE") %>%
    tibble::as_tibble() %>%
    dplyr::filter(SampleID %in% data_features$SampleID) %>%
    dplyr::filter(feature %in% names(dplyr::select(data_features, contains("OID")))) %>%
    dplyr::select(SampleID, feature, LOD) %>%
    unique()
  
  # label ====
  LABEL <- basename(i) %>%
    gsub(pattern = "\\.rda$", replacement = "")
  cat("Processing: ", LABEL, "\n")
  
  # make analysis df ====
  df <- data_features %>%
    dplyr::left_join(data_meta_samples, by = "SampleID") %>%
    dplyr::mutate(
      Smoke_Stat = factor(Smoke_Stat, levels = c("Never", "Former", "Smoker", "Unknown")),
      Pa_Index = factor(Pa_Index, levels = c("Inactive", "Moderately inactive", "Moderately active", "Active", "Missing")),
      L_School = factor(L_School, levels = c("None", "Primary", "Secondary", "Technical/professional", "Longer education", "Not specified"))
    ) 
  
  # make analysis list ====
  site_proximal <- c("C180", "C181", "C182", "C183", "C184", "C185")
  site_distal <- c("C186", "C187")
  site_exclude <- c("C188", "C189")
  
  ## lists ====
  list <- list(
    # overall
    overall_combined_0_heterogeneity = {
      temp <- df %>%
        dplyr::filter(!Siteclrt %in% site_exclude)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp) %>%
        dplyr::mutate(
          heterogeneity_site_temp = case_when(
            Siteclrt %in% site_proximal ~ "proximal",
            Siteclrt %in% site_distal ~ "distal",
            TRUE ~ NA_character_)) %>%
        dplyr::group_by(Match_Caseset) %>%
        dplyr::mutate(
          heterogeneity_site = first(na.omit(heterogeneity_site_temp))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-heterogeneity_site_temp)
    },
    
    # overall 2y
    overall_combined_2_heterogeneity = {
      temp <- df %>%
        dplyr::filter(!Siteclrt %in% site_exclude) %>%
        dplyr::mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        dplyr::filter(Length_Bld_Modified > 2)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp) %>%
        dplyr::mutate(
          heterogeneity_site_temp = case_when(
            Siteclrt %in% site_proximal ~ "proximal",
            Siteclrt %in% site_distal ~ "distal",
            TRUE ~ NA_character_)) %>%
        dplyr::group_by(Match_Caseset) %>%
        dplyr::mutate(
          heterogeneity_site = first(na.omit(heterogeneity_site_temp))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-heterogeneity_site_temp)
    },
    
    # overall 5y
    overall_combined_5_heterogeneity = {
      temp <- df %>%
        dplyr::filter(!Siteclrt %in% site_exclude) %>%
        dplyr::mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        dplyr::filter(Length_Bld_Modified > 5)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp) %>%
        dplyr::mutate(
          heterogeneity_site_temp = case_when(
            Siteclrt %in% site_proximal ~ "proximal",
            Siteclrt %in% site_distal ~ "distal",
            TRUE ~ NA_character_)) %>%
        dplyr::group_by(Match_Caseset) %>%
        dplyr::mutate(
          heterogeneity_site = first(na.omit(heterogeneity_site_temp))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-heterogeneity_site_temp)
    }
  )
  
  # analysis ====
  data_results <- analysis_heterogeneity(data_list = list, 
                                         models = models, 
                                         models_heterogeneity = models_heterogeneity, 
                                         exposure = "OID", 
                                         outcome_var = "Cncr_Caco_Clrt",
                                         heterogeneity_col = "Sex", 
                                         strata_var = "Match_Caseset", 
                                         base_covariates = "Age_Recr",
                                         covariates = "Bmi_C + Alc_Re + Smoke_Stat + Pa_Index + L_School")
  write.table(data_results, paste0("analysis/002_EPIC-analysis/002_heterogeneity/", VAR[2], "_", VAR[3], "/heterogeneity-site_", LABEL, ".txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

  data_results <- analysis_heterogeneity(data_list = list, 
                                         models = models, 
                                         models_heterogeneity = models_heterogeneity, 
                                         exposure = "OID", 
                                         outcome_var = "Cncr_Caco_Clrt",
                                         heterogeneity_col = "heterogeneity_site", 
                                         strata_var = "Match_Caseset", 
                                         base_covariates = "Age_Recr",
                                         covariates = "Bmi_C + Alc_Re + Smoke_Stat + Pa_Index + L_School")
  write.table(data_results, paste0("analysis/002_EPIC-analysis/002_heterogeneity/", VAR[2], "_", VAR[3], "/heterogeneity-sex_", LABEL, ".txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

# Parallel setup ====
cl <- parallel::makeCluster(spec = 8)  # Create the cluster
# Step 1: Source and load packages on each worker
parallel::clusterEvalQ(cl, {
  source("src/000_source.R")
  source("src/000_functions.R")
})
# Step 2: Export main function
parallel::clusterExport(cl, varlist = c("process_file"), envir = environment())
# Step 3: Run parallel processing
results <- parallel::parLapply(cl, LIST_FILES, process_file)

# Step 4: Clean up
stopCluster(cl)
