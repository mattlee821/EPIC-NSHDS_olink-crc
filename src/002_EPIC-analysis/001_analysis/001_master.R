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
  
  list <- list(
    # overall
    overall_combined_0 = df,
    overall_male_0 = df %>% dplyr::filter(Sex == 1),
    overall_female_0 = df %>% dplyr::filter(Sex == 2),
    # proximal
    proximal_combined_0 = {
      temp <- df %>% dplyr::filter(Siteclrt %in% site_proximal)
      temp <- unique(temp$Match_Caseset)
      df %>% dplyr::filter(Match_Caseset %in% temp)
    },
    proximal_male_0 = {
      temp <- df %>% dplyr::filter(Siteclrt %in% site_proximal)
      temp <- unique(temp$Match_Caseset)
      df %>% 
        dplyr::filter(Match_Caseset %in% temp) %>%
        dplyr::filter(Sex == 1)
    },
    proximal_female_0 = {
      temp <- df %>% dplyr::filter(Siteclrt %in% site_proximal)
      temp <- unique(temp$Match_Caseset)
      df %>% 
        dplyr::filter(Match_Caseset %in% temp) %>%
        dplyr::filter(Sex == 2)
    },
    # distal
    distal_combined_0 = {
      temp <- df %>% dplyr::filter(Siteclrt %in% site_distal)
      temp <- unique(temp$Match_Caseset)
      df %>% dplyr::filter(Match_Caseset %in% temp)
    },
    distal_male_0 = {
      temp <- df %>% dplyr::filter(Siteclrt %in% site_distal)
      temp <- unique(temp$Match_Caseset)
      df %>% 
        dplyr::filter(Match_Caseset %in% temp) %>%
        dplyr::filter(Sex == 1)
    },
    distal_female_0 = {
      temp <- df %>% dplyr::filter(Siteclrt %in% site_distal)
      temp <- unique(temp$Match_Caseset)
      df %>% 
        dplyr::filter(Match_Caseset %in% temp) %>%
        dplyr::filter(Sex == 2)
    }
  )
  
  # make analysis list excluding 2y follow up ====
  site_proximal <- c("C180", "C181", "C182", "C183", "C184", "C185")
  site_distal <- c("C186", "C187")
  exclusion <- 2
  list_2y <- list(
    # overall
    overall_combined_2 = {
      temp <- df %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    overall_male_2 = {
      temp <- df %>% dplyr::filter(Sex == 1) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    overall_female_2 = {
      temp <- df %>% dplyr::filter(Sex == 2) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    # proximal
    proximal_combined_2 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_proximal) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    proximal_male_2 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_proximal) %>%
        dplyr::filter(Sex == 1) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    proximal_female_2 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_proximal) %>%
        dplyr::filter(Sex == 2) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    # distal 
    distal_combined_2 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_distal) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    distal_male_2 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_distal) %>%
        dplyr::filter(Sex == 1) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    distal_female_2 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_distal) %>%
        dplyr::filter(Sex == 2) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    }
  )
  
  # make analysis list excluding 5y follow up ====
  site_proximal <- c("C180", "C181", "C182", "C183", "C184", "C185")
  site_distal <- c("C186", "C187")
  exclusion <- 5
  list_5y <- list(
    # overall
    overall_combined_5 = {
      temp <- df %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    overall_male_5 = {
      temp <- df %>% dplyr::filter(Sex == 1) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    overall_female_5 = {
      temp <- df %>% dplyr::filter(Sex == 2) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    # proximal
    proximal_combined_5 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_proximal) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    proximal_male_5 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_proximal) %>%
        dplyr::filter(Sex == 1) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    proximal_female_5 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_proximal) %>%
        dplyr::filter(Sex == 2) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    # distal 
    distal_combined_5 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_distal) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    distal_male_5 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_distal) %>%
        dplyr::filter(Sex == 1) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    },
    distal_female_5 = {
      temp <- df %>%
        dplyr::filter(Siteclrt %in% site_distal) %>%
        dplyr::filter(Sex == 2) %>%
        mutate(Length_Bld_Modified = Length_Bld / 365) %>%
        filter(Length_Bld_Modified > exclusion)
      temp <- temp %>%
        filter(Cncr_Caco_Clrt == 1)
      temp <- temp %>%
        pull(Match_Caseset) %>%
        unique()
      df %>% filter(Match_Caseset %in% temp)
    }
  )
  
  # analysis ====
  list <- c(list, list_2y, list_5y)
  data_results <- analysis(data_list = list, models = models, covariates = "Bmi_C + Alc_Re + Smoke_Stat + Pa_Index + L_School")
  
  write.table(data_results, paste0("analysis/002_EPIC-analysis/001_analysis/", VAR_data, "/", LABEL, ".txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
