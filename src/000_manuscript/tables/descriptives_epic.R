rm(list=ls())
set.seed(821)

# descriptive table 

# environment ====
source("src/000_source.R")
rm(list=ls())

# data_raw ====
# data_samples
data_samples <- haven::read_sas(data_file = "data/raw/EPIC_phenotype/clrt_caco.sas7bdat")
temp_mapping <- data_samples %>%
  dplyr::select(Idepic_Bio) %>%
  tidyr::extract(Idepic_Bio, into = c("prefix", "SampleID"), regex = "^([^_]+)_(.+)$", remove = FALSE) %>%
  dplyr::mutate(SampleID = gsub("_", "", SampleID)) %>%
  dplyr::select(-prefix)
# data olink
data_immuneonc <- OlinkAnalyze::read_NPX(filename = "data/raw/EPIC_olink-immuneonc/20211440_Murphy_NPX.xlsx") %>%
  tibble::as_tibble()
data_explore <- OlinkAnalyze::read_NPX(filename = "data/raw/EPIC_olink-explore/WL-3904_NPX_2024-12-09.csv") %>%
  tibble::as_tibble()
# formatting
data_immuneonc <- data_immuneonc %>%
  dplyr::mutate(Sample_Type = ifelse(grepl("CONTROL", SampleID), "CONTROL", "SAMPLE"))
data_explore <- data_explore %>%
  dplyr::mutate(Sample_Type = ifelse(SampleID == "Empty well", "CONTROL", Sample_Type)) %>%
  dplyr::mutate(
    SampleID = if_else(
      Sample_Type == "SAMPLE",
      stringr::str_pad(as.character(as.numeric(SampleID)), width = 5, pad = "0"),
      SampleID
    )
  ) %>%
  dplyr::left_join(temp_mapping, by = "SampleID") %>%
  dplyr::mutate(SampleID = ifelse(Sample_Type == "SAMPLE", Idepic_Bio, SampleID)) %>%
  dplyr::select(-Idepic_Bio)
# pheno join
data_immuneonc <- data_immuneonc %>%
  dplyr::left_join(data_samples, by = c("SampleID" = "Idepic_Bio"))
data_explore <- data_explore %>%
  dplyr::left_join(data_samples, by = c("SampleID" = "Idepic_Bio"))
# VARS
VAR_N_immuneonc_raw <- data_immuneonc %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(SampleID) %>% unique() %>% length()
VAR_N_explore_raw <- data_explore %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(SampleID) %>% unique() %>% length()
VAR_caseset_immuneonc_raw <- data_immuneonc %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(Match_Caseset) %>% unique() %>% length()
VAR_caseset_explore_raw <- data_explore %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(Match_Caseset) %>% unique() %>% length()
VAR_feature_immuneonc_raw <- data_immuneonc %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(OlinkID) %>% unique() %>% length()
VAR_feature_explore_raw <- data_explore %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(OlinkID) %>% unique() %>% length()

# data_processed1 ====
data_immuneonc <- data.table::fread("data/processed/EPIC_olink-immuneonc.txt")
data_explore <- data.table::fread("data/processed/EPIC_olink-explore.txt")
# VARS
VAR_N_immuneonc_processed1 <- data_immuneonc %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(SampleID) %>% unique() %>% length()
VAR_N_explore_processed1 <- data_explore %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(SampleID) %>% unique() %>% length()
VAR_caseset_immuneonc_processed1 <- data_immuneonc %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(Match_Caseset) %>% unique() %>% length()
VAR_caseset_explore_processed1 <- data_explore %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(Match_Caseset) %>% unique() %>% length()
VAR_feature_immuneonc_processed1 <- data_immuneonc %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(OlinkID) %>% unique() %>% length()
VAR_feature_explore_processed1 <- data_explore %>% dplyr::filter(Sample_Type == "SAMPLE") %>% dplyr::pull(OlinkID) %>% unique() %>% length()

# data_processed2 ====
load("analysis/001_data-processing/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.rda")
data_immuneonc <- df_processed
load("analysis/001_data-processing/EPIC_olink-explore_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.rda")
data_explore <- df_processed
# VARS
VAR_N_immuneonc_processed2 <- data_immuneonc %>% dplyr::pull(SampleID) %>% unique() %>% length()
VAR_N_explore_processed2 <- data_explore %>% dplyr::pull(SampleID) %>% unique() %>% length()
VAR_feature_immuneonc_processed2 <- sum(grepl("OID", colnames(data_immuneonc)))
VAR_feature_explore_processed2 <- sum(grepl("OID", colnames(data_explore)))

# table ====
table_N <- data.frame(
  data = c("immune-oncology", "explore", "immune-oncology", "explore", "immune-oncology", "explore"),
  step = c("raw", "raw", "processed 1", "processed 1", "processed 2", "processed 2"),
  samples = c(VAR_N_immuneonc_raw, VAR_N_explore_raw, VAR_N_immuneonc_processed1, VAR_N_explore_processed1, VAR_N_immuneonc_processed2, VAR_N_explore_processed2),
  case_sets = c(VAR_caseset_immuneonc_raw, VAR_caseset_explore_raw, VAR_caseset_immuneonc_processed1, VAR_caseset_explore_processed1, VAR_caseset_immuneonc_processed2, VAR_caseset_explore_processed2),
  features = c(VAR_feature_immuneonc_raw, VAR_feature_explore_raw, VAR_feature_immuneonc_processed1, VAR_feature_explore_processed1, VAR_feature_immuneonc_processed2, VAR_feature_explore_processed2)
)
write.table(table_N, "manuscript/tables/N-processing_epic.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# descriptives ====
load("analysis/001_data-processing/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.rda")
id_immuneonc_samples <- df_processed %>%
  pull(SampleID)
id_immuneonc_features <- colnames(df_processed)[grepl("OID", colnames(df_processed))]
load("analysis/001_data-processing/EPIC_olink-explore_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.rda")
id_explore_samples <- df_processed %>%
  pull(SampleID)
id_explore_features <- colnames(df_processed)[grepl("OID", colnames(df_processed))]

data_immuneonc <- data.table::fread("data/processed/EPIC_olink-immuneonc.txt") %>%
  filter(SampleID %in% id_immuneonc_samples) %>%
  filter(OlinkID %in% id_immuneonc_features) %>%
  tibble:::as_tibble()
data_explore <- data.table::fread("data/processed/EPIC_olink-explore.txt") %>%
  filter(SampleID %in% id_explore_samples) %>%
  filter(OlinkID %in% id_explore_features) %>%
  tibble:::as_tibble()

# functions ====
compute_descriptives <- function(data) {
  data %>%
    summarise(
      N = n(),
      age_median = median(Age_Recr, na.rm = TRUE),
      age_sd = sd(Age_Recr, na.rm = TRUE),
      followup_median = median(Agexit - Age_Recr, na.rm = TRUE),
      followup_sd = sd(Agexit - Age_Recr, na.rm = TRUE),
      followup_iqr = IQR(Agexit - Age_Recr, na.rm = TRUE),
      BMI_median = median(Bmi_C, na.rm = TRUE),
      BMI_sd = sd(Bmi_C, na.rm = TRUE),
      alcohol_median = median(Alc_Re, na.rm = TRUE),
      alcohol_sd = sd(Alc_Re, na.rm = TRUE)
    )
}
compute_categorical <- function(data) {
  list(
    smoking = table(data$Smoke_Stat),
    physical_activity = table(data$Pa_Index),
    education = table(data$L_School)
  )
}
desc_to_long <- function(df, name) {
  df %>%
    dplyr::ungroup() %>%                               # remove grouping metadata
    dplyr::select(where(~ is.numeric(.x))) %>%         # keep only numeric columns
    tidyr::pivot_longer(everything(), names_to = "statistic", values_to = name)
}
cat_to_long <- function(cat_list, name) {
  dplyr::bind_rows(
    purrr::imap(cat_list, function(tbl, var_name) {
      as.data.frame(tbl) %>%
        dplyr::rename(level = 1, value = Freq) %>%
        dplyr::mutate(statistic = paste0(var_name, ": ", level)) %>%
        dplyr::select(statistic, value)
    })
  ) %>%
    dplyr::group_by(statistic) %>%
    dplyr::summarise(!!name := value, .groups = "drop")
}
# immuneonc ====
site_proximal <- c("C180", "C181", "C182", "C183", "C184", "C185")
site_distal <- c("C186", "C187")

# Recode Cncr_Caco_Clrt to "Case"/"Control"
data_phenofile <- data_immuneonc %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  mutate(Cncr_Caco_Clrt = dplyr::recode(as.character(Cncr_Caco_Clrt),
                                        `0` = "Control", `1` = "Case"),
         Sex = dplyr::recode(as.character(Sex),
                             `1` = "Male", `2` = "Female"))

# For all samples
desc_all <- compute_descriptives(data_phenofile)
cat_all <- compute_categorical(data_phenofile)

# By sex
desc_by_sex <- data_phenofile %>%
  group_by(Sex) %>%
  group_modify(~ compute_descriptives(.x))

cat_by_sex <- data_phenofile %>%
  group_by(Sex) %>%
  group_split() %>%
  setNames(paste0("Sex_", sapply(., function(x) unique(x$Sex)))) %>%
  lapply(compute_categorical)

# By case/control (now with Case/Control labels)
desc_by_case <- data_phenofile %>%
  group_by(Cncr_Caco_Clrt) %>%
  group_modify(~ compute_descriptives(.x))

cat_by_case <- data_phenofile %>%
  group_by(Cncr_Caco_Clrt) %>%
  group_split() %>%
  setNames(paste0("Case_", sapply(., function(x) unique(x$Cncr_Caco_Clrt)))) %>%
  lapply(compute_categorical)

# By sex and case/control
desc_by_sex_case <- data_phenofile %>%
  group_by(Sex, Cncr_Caco_Clrt) %>%
  group_modify(~ compute_descriptives(.x))

cat_by_sex_case <- data_phenofile %>%
  group_by(Sex, Cncr_Caco_Clrt) %>%
  group_split() %>%
  setNames(nm = data_phenofile %>%
             group_by(Sex, Cncr_Caco_Clrt) %>%
             group_keys() %>%
             mutate(name = paste0("Sex_", Sex, "_Case_", Cncr_Caco_Clrt)) %>%
             pull(name)) %>%
  lapply(compute_categorical)

# Matched proximal
proximal_cases <- data_phenofile %>%
  filter(Cncr_Caco_Clrt == "Case", Siteclrt %in% site_proximal)

proximal_controls <- data_phenofile %>%
  filter(Cncr_Caco_Clrt == "Control", Match_Caseset %in% proximal_cases$Match_Caseset)

desc_prox_cases <- compute_descriptives(proximal_cases)
desc_prox_controls <- compute_descriptives(proximal_controls)
cat_prox_cases <- compute_categorical(proximal_cases)
cat_prox_controls <- compute_categorical(proximal_controls)

# Matched distal
distal_cases <- data_phenofile %>%
  filter(Cncr_Caco_Clrt == "Case", Siteclrt %in% site_distal)

distal_controls <- data_phenofile %>%
  filter(Cncr_Caco_Clrt == "Control", Match_Caseset %in% distal_cases$Match_Caseset)

desc_dist_cases <- compute_descriptives(distal_cases)
desc_dist_controls <- compute_descriptives(distal_controls)
cat_dist_cases <- compute_categorical(distal_cases)
cat_dist_controls <- compute_categorical(distal_controls)

# Create long dataframes
desc_tables <- list(
  All = desc_all,
  Sex_Male = dplyr::filter(desc_by_sex, Sex == "Male"),
  Sex_Female = dplyr::filter(desc_by_sex, Sex == "Female"),
  Case = dplyr::filter(desc_by_case, Cncr_Caco_Clrt == "Case"),
  Control = dplyr::filter(desc_by_case, Cncr_Caco_Clrt == "Control"),
  Sex_Male_Case = dplyr::filter(desc_by_sex_case, Sex == "Male", Cncr_Caco_Clrt == "Case"),
  Sex_Female_Case = dplyr::filter(desc_by_sex_case, Sex == "Female", Cncr_Caco_Clrt == "Case"),
  Sex_Male_Control = dplyr::filter(desc_by_sex_case, Sex == "Male", Cncr_Caco_Clrt == "Control"),
  Sex_Female_Control = dplyr::filter(desc_by_sex_case, Sex == "Female", Cncr_Caco_Clrt == "Control"),
  Proximal_Case = desc_prox_cases,
  Proximal_Control = desc_prox_controls,
  Distal_Case = desc_dist_cases,
  Distal_Control = desc_dist_controls
)

desc_long <- purrr::imap(desc_tables, desc_to_long) %>%
  purrr::reduce(full_join, by = "statistic")

cat_tables <- list(
  All = cat_all,
  Sex_Male = cat_by_sex$Sex_Male,
  Sex_Female = cat_by_sex$Sex_Female,
  Case = cat_by_case$Case_Case,
  Control = cat_by_case$Case_Control,
  Sex_Male_Case = cat_by_sex_case$Sex_Male_Case_Case,
  Sex_Female_Case = cat_by_sex_case$Sex_Female_Case_Case,
  Sex_Male_Control = cat_by_sex_case$Sex_Male_Case_Control,
  Sex_Female_Control = cat_by_sex_case$Sex_Female_Case_Control,
  Proximal_Case = cat_prox_cases,
  Proximal_Control = cat_prox_controls,
  Distal_Case = cat_dist_cases,
  Distal_Control = cat_dist_controls
)

cat_long <- purrr::imap(cat_tables, cat_to_long) %>%
  purrr::reduce(full_join, by = "statistic")

# combine
desired_order <- c(
  "N",
  "age_median", 
  "age_sd", 
  "followup_median", 
  "followup_sd", 
  "followup_iqr", 
  "BMI_median", 
  "BMI_sd", 
  "alcohol_median", 
  "alcohol_sd",
  "smoking: Never",
  "smoking: Former",
  "smoking: Smoker",
  "smoking: Unknown",
  "physical_activity: Inactive",
  "physical_activity: Moderately inactive",
  "physical_activity: Moderately active",
  "physical_activity: Active",
  "physical_activity: Missing",
  "education: None",
  "education: Primary",
  "education: Secondary",
  "education: Technical/professional",
  "education: Longer education",
  "education: Not specified"
)

summary_table <- dplyr::bind_rows(desc_long, cat_long) %>%
  dplyr::mutate(statistic = factor(statistic, levels = desired_order)) %>%
  dplyr::arrange(statistic) %>%
  dplyr::select(statistic, All, Case, Control, Sex_Male, Sex_Male_Case, Sex_Male_Control, Sex_Female, Sex_Female_Case, Sex_Female_Control, Proximal_Case, Proximal_Control, Distal_Case, Distal_Control) %>%
  dplyr::rename_with(~ gsub("_?Sex_?|_?$", "", .x)) %>%
  dplyr::mutate(
    statistic = as.character(statistic),
    dplyr::across(
      where(is.numeric),
      ~ case_when(
        statistic %in% c(
          "N",
          "smoking: never", "smoking: former", "smoking: smoker", "smoking: unknown",
          "physical_activity: Inactive", "physical_activity: Moderately inactive",
          "physical_activity: Moderately active", "physical_activity: Active", "physical_activity: Missing",
          "education: None", "education: Primary", "education: Secondary",
          "education: Technical/professional", "education: Longer education", "education: Not specified"
        ) ~ as.numeric(round(.x)),  # ensures no decimal
        TRUE ~ round(.x, 2)
      )
    )
  )

# write
write.table(summary_table, "manuscript/tables/descriptives_epic-immuneonc.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# explore ====
site_proximal <- c("C180", "C181", "C182", "C183", "C184", "C185")
site_distal <- c("C186", "C187")

# Recode Cncr_Caco_Clrt to "Case"/"Control"
data_phenofile <- data_explore %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  mutate(Cncr_Caco_Clrt = dplyr::recode(as.character(Cncr_Caco_Clrt),
                                        `0` = "Control", `1` = "Case"))
# For all samples
desc_all <- compute_descriptives(data_phenofile)
cat_all <- compute_categorical(data_phenofile)

# By sex
desc_by_sex <- data_phenofile %>%
  group_by(Sex) %>%
  group_modify(~ compute_descriptives(.x))

cat_by_sex <- data_phenofile %>%
  group_by(Sex) %>%
  group_split() %>%
  setNames(paste0("Sex_", sapply(., function(x) unique(x$Sex)))) %>%
  lapply(compute_categorical)

# By case/control (now with Case/Control labels)
desc_by_case <- data_phenofile %>%
  group_by(Cncr_Caco_Clrt) %>%
  group_modify(~ compute_descriptives(.x))

cat_by_case <- data_phenofile %>%
  group_by(Cncr_Caco_Clrt) %>%
  group_split() %>%
  setNames(paste0("Case_", sapply(., function(x) unique(x$Cncr_Caco_Clrt)))) %>%
  lapply(compute_categorical)

# By sex and case/control
desc_by_sex_case <- data_phenofile %>%
  group_by(Sex, Cncr_Caco_Clrt) %>%
  group_modify(~ compute_descriptives(.x))

cat_by_sex_case <- data_phenofile %>%
  group_by(Sex, Cncr_Caco_Clrt) %>%
  group_split() %>%
  setNames(nm = data_phenofile %>%
             group_by(Sex, Cncr_Caco_Clrt) %>%
             group_keys() %>%
             mutate(name = paste0("Sex_", Sex, "_Case_", Cncr_Caco_Clrt)) %>%
             pull(name)) %>%
  lapply(compute_categorical)

# Matched proximal
proximal_cases <- data_phenofile %>%
  filter(Cncr_Caco_Clrt == "Case", Siteclrt %in% site_proximal)

proximal_controls <- data_phenofile %>%
  filter(Cncr_Caco_Clrt == "Control", Match_Caseset %in% proximal_cases$Match_Caseset)

desc_prox_cases <- compute_descriptives(proximal_cases)
desc_prox_controls <- compute_descriptives(proximal_controls)
cat_prox_cases <- compute_categorical(proximal_cases)
cat_prox_controls <- compute_categorical(proximal_controls)

# Matched distal
distal_cases <- data_phenofile %>%
  filter(Cncr_Caco_Clrt == "Case", Siteclrt %in% site_distal)

distal_controls <- data_phenofile %>%
  filter(Cncr_Caco_Clrt == "Control", Match_Caseset %in% distal_cases$Match_Caseset)

desc_dist_cases <- compute_descriptives(distal_cases)
desc_dist_controls <- compute_descriptives(distal_controls)
cat_dist_cases <- compute_categorical(distal_cases)
cat_dist_controls <- compute_categorical(distal_controls)

# Create long dataframes
desc_tables <- list(
  All = desc_all,
  Sex_Male = dplyr::filter(desc_by_sex, Sex == "Male"),
  Sex_Female = dplyr::filter(desc_by_sex, Sex == "Female"),
  Case = dplyr::filter(desc_by_case, Cncr_Caco_Clrt == "Case"),
  Control = dplyr::filter(desc_by_case, Cncr_Caco_Clrt == "Control"),
  Sex_Male_Case = dplyr::filter(desc_by_sex_case, Sex == "Male", Cncr_Caco_Clrt == "Case"),
  Sex_Female_Case = dplyr::filter(desc_by_sex_case, Sex == "Female", Cncr_Caco_Clrt == "Case"),
  Sex_Male_Control = dplyr::filter(desc_by_sex_case, Sex == "Male", Cncr_Caco_Clrt == "Control"),
  Sex_Female_Control = dplyr::filter(desc_by_sex_case, Sex == "Female", Cncr_Caco_Clrt == "Control"),
  Proximal_Case = desc_prox_cases,
  Proximal_Control = desc_prox_controls,
  Distal_Case = desc_dist_cases,
  Distal_Control = desc_dist_controls
)

desc_long <- purrr::imap(desc_tables, desc_to_long) %>%
  purrr::reduce(full_join, by = "statistic")

cat_tables <- list(
  All = cat_all,
  Sex_Male = cat_by_sex$Sex_Male,
  Sex_Female = cat_by_sex$Sex_Female,
  Case = cat_by_case$Case_Case,
  Control = cat_by_case$Case_Control,
  Sex_Male_Case = cat_by_sex_case$Sex_Male_Case_Case,
  Sex_Female_Case = cat_by_sex_case$Sex_Female_Case_Case,
  Sex_Male_Control = cat_by_sex_case$Sex_Male_Case_Control,
  Sex_Female_Control = cat_by_sex_case$Sex_Female_Case_Control,
  Proximal_Case = cat_prox_cases,
  Proximal_Control = cat_prox_controls,
  Distal_Case = cat_dist_cases,
  Distal_Control = cat_dist_controls
)

cat_long <- purrr::imap(cat_tables, cat_to_long) %>%
  purrr::reduce(full_join, by = "statistic")

# combine
desired_order <- c(
  "N",
  "age_median", 
  "age_sd", 
  "followup_median", 
  "followup_sd", 
  "followup_iqr", 
  "BMI_median", 
  "BMI_sd", 
  "alcohol_median", 
  "alcohol_sd",
  "smoking: Never",
  "smoking: Former",
  "smoking: Smoker",
  "smoking: Unknown",
  "physical_activity: Inactive",
  "physical_activity: Moderately inactive",
  "physical_activity: Moderately active",
  "physical_activity: Active",
  "physical_activity: Missing",
  "education: None",
  "education: Primary",
  "education: Secondary",
  "education: Technical/professional",
  "education: Longer education",
  "education: Not specified"
)

summary_table <- dplyr::bind_rows(desc_long, cat_long) %>%
  dplyr::mutate(statistic = factor(statistic, levels = desired_order)) %>%
  dplyr::arrange(statistic) %>%
  dplyr::select(statistic, All, Case, Control, Sex_Male, Sex_Male_Case, Sex_Male_Control, Sex_Female, Sex_Female_Case, Sex_Female_Control, Proximal_Case, Proximal_Control, Distal_Case, Distal_Control) %>%
  dplyr::rename_with(~ gsub("_?Sex_?|_?$", "", .x)) %>%
  dplyr::mutate(
    statistic = as.character(statistic),
    dplyr::across(
      where(is.numeric),
      ~ case_when(
        statistic %in% c(
          "N",
          "smoking: never", "smoking: former", "smoking: smoker", "smoking: unknown",
          "physical_activity: Inactive", "physical_activity: Moderately inactive",
          "physical_activity: Moderately active", "physical_activity: Active", "physical_activity: Missing",
          "education: None", "education: Primary", "education: Secondary",
          "education: Technical/professional", "education: Longer education", "education: Not specified"
        ) ~ as.numeric(round(.x)),  # ensures no decimal
        TRUE ~ round(.x, 2)
      )
    )
  )

# write
write.table(summary_table, "manuscript/tables/descriptives_epic-explore.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
