rm(list=ls())
set.seed(821)

#environment ====
library(dplyr)
library(survival)
source("src/000_functions.R")

# files ====
LIST_FILES <- list.files("analysis/001_data-processing/", pattern = "data-features_", all.files = T, full.names = T, recursive = T)
VAR_data <- "EPIC_olink-immuneonc_data-processed_LOD-filter"
LIST_FILES <- LIST_FILES[grepl(VAR_data, LIST_FILES)]
df_meta_samples <- data.table::fread(paste0("data/processed/EPIC_olink-immuneonc_sample-metadata.txt"))

i = LIST_FILES[1]

for (i in LIST_FILES){
  # data ====
  data_meta_samples <- df_meta_samples
  LABEL <- i
  load(LABEL)
  LABEL <- basename(LABEL)
  LABEL <- gsub(".rda", "", LABEL)
  data_features <- df_processed
  cat("Processing: ", LABEL, "\n")
  
  # filter data_meta_samples ====
  data_meta_samples <- data_meta_samples %>%
    dplyr::filter(Idepic_Bio %in% data_features$Idepic_Bio)
  
  # make analysis df ====
  df <- data_features %>%
    dplyr::left_join(data_meta_samples, by = "Idepic_Bio")
  
  # format df ====
  df <- df %>%
    # make age groups
    dplyr::mutate(
      age_group = as.integer(cut(Age_Recr,
                                 breaks = seq(0, max(Age_Recr) + 5, by = 5),
                                 right = FALSE,
                                 labels = seq(0, max(Age_Recr), by = 5))),
      
      # make day of blood draw
      day_of_year = data.table::yday(as.Date(D_Bld_Coll, format = "%Y-%m-%d")), # Convert date of blood collection to day of year
      ## Trigonometric transformations
      ## Regarding adjustment for month/season of blood draw, in analyses of circulating vitamin D we have had success with spherical adjustment for day-of-year of blood draw. 
      ## By using a (very) truncated Fourier transformation, I think we will be able to get away with four parameters for this adjustment rather than the 11 parameters needed to adjust for month. 
      ## This also has the added bonus that the adjustment is smooth, with no artificial discontinuities from season to season or year to year.
      day_of_year_blood_draw_sin_2_pi = sin(2 * pi * day_of_year / 365),
      day_of_year_blood_draw_cos_2_pi = cos(2 * pi * day_of_year / 365),
      day_of_year_blood_draw_sin_4_pi = sin(4 * pi * day_of_year / 365),
      day_of_year_blood_draw_cos_4_pi = cos(4 * pi * day_of_year / 365),
      
      # covariate formatting
      # sex: 1 = male; 2 = female
      # smoking status: 1 = never; 2 = former; 3 = smoker; 4 = unknown
      Smoke_Stat = as.factor(Smoke_Stat),
      # physical activity index: 1 = inactive; 2 = moderately inactive; 3 = moderately active; 4 = active; 5 = missing
      Pa_Index = as.factor(Pa_Index),
      # level of schooling: 0 = none; 1 = primary; 2 = technical/professional; 3 = secondary; 4 = longer; 5 = not specified
      L_School = as.factor(L_School)
      
      # factor ordering
      # Smoke_Stat = factor(Smoke_Stat, levels = c("Never", "Former", "Smoker")),
      # Pa_Index = factor(Pa_Index, levels = c("Inactive", "Moderately inactive", "Moderately active", "Active", "Missing")),
      # L_School = factor(L_School, levels = c("None", "Primary school completed", "Secondary school", "Technical/professional school", "Longer education (incl. University deg.)", "Not specified"))
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
  data_results <- analysis_heterogeneity(data_list = list, models = models, models_heterogeneity = models_heterogeneity, heterogeneity_col = "Sex", covariates = "Bmi_C + Alc_Re + Smoke_Stat + Pa_Index + L_School")
  write.table(data_results, paste0("analysis/002_EPIC-analysis/002_heterogeneity/", VAR_data, "/heterogeneity-site_", LABEL, ".txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  data_results <- analysis_heterogeneity(data_list = list, models = models, models_heterogeneity = models_heterogeneity, heterogeneity_col = "heterogeneity_site", covariates = "Bmi_C + Alc_Re + Smoke_Stat + Pa_Index + L_School")
  write.table(data_results, paste0("analysis/002_EPIC-analysis/002_heterogeneity/", VAR_data, "/heterogeneity-sex_", LABEL, ".txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
