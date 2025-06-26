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
  data_results <- analysis(data_list = list, models = models, 
                           exposure = "OID", 
                           outcome_var = "Cncr_Caco_Clrt",
                           strata_var = "Match_Caseset", 
                           base_covariates = "Age_Recr",
                           covariates = "Bmi_C + Alc_Re + Smoke_Stat + Pa_Index + L_School")
  
  write.table(data_results, paste0("analysis/002_EPIC-analysis/001_analysis/", VAR[2], "_", VAR[3], "/", LABEL, ".txt"), 
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
