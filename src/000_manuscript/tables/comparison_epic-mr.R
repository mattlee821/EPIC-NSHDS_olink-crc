rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(openxlsx)
library(data.table)
library(stringr)

# EPIC/MR comparison ====
## data epic ====
data <- data.table::fread("analysis/002_coxph/003_format/001_analysis/cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP/formatted-analysis.txt")
data <- data %>%
  filter(cancer != "earlyonset",
         Organism == "Human") %>%
  mutate(cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
         sex = factor(sex, levels = c("combined", "female", "male")),
         followup = factor(followup, levels = c(0, 2, 5)),
         raw_complete = factor(raw_complete, levels = c("raw", "complete"))) %>%
  arrange(cancer, sex, followup, raw_complete, data_id, Target, SeqId, model_type) %>%
  group_by(cancer, sex, raw_complete, followup, model_type, Organism) %>%
  mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  ungroup() %>%
  filter(
    followup == 0,
    raw_complete == "raw",
    model_type == "model_2_extra",
    FDR_BH < 0.05
  ) %>%
  select(SeqId, UniProt, Target, TargetFullName, cancer, sex, coef, se_robust, pval, FDR_BH)
ID_seqid <- data$SeqId %>%
  unique() %>%
  str_replace_all("-", "_")
ID_uniprot <- unique(data$UniProt) 

# data_mr forward ====
LIST_FILES <- list.files("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/001_exposure-proteins_outcome-cancer/MR/", pattern = ".txt", all.files = T, full.names = T)
LIST_FILES <- LIST_FILES[str_detect(LIST_FILES, str_c(c(ID_seqid, ID_uniprot), collapse = "|"))]
data_mr <- lapply(LIST_FILES, data.table::fread, header = T, sep = "\t")
data_mr <- bind_rows(data_mr)
## remove all coloc.abf 
data_mr <- data_mr %>%
  filter(!(exposure_finemap == "coloc.abf")) 
## add in all singlesnp MR for coloc.abf
LIST_FILES <- list.files("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/001_exposure-proteins_outcome-cancer/singlesnp/", pattern = ".txt", all.files = T, full.names = T)
LIST_FILES <- LIST_FILES[str_detect(LIST_FILES, str_c(c(ID_seqid, ID_uniprot), collapse = "|"))]
data_mrsinglesnp <- lapply(LIST_FILES, data.table::fread, header = T, sep = "\t")
data_mrsinglesnp <- bind_rows(data_mrsinglesnp)
data_mrsinglesnp <- data_mrsinglesnp %>%
  filter(exposure_finemap == "coloc.abf") %>%
  mutate(method = "Wald ratio")
data_mr <- data_mr %>%
  bind_rows(data_mrsinglesnp)
## format for seqid/uniprot matching
data_mr_forward <- data_mr %>%
  mutate(
    SeqId = case_when(
      exposure_study == "ferkingstad_2021_PMID34857953" ~ sub("^([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", exposure),
      exposure_study == "zhang_2022_PMID35501419" ~ sub("^SeqId_(.*)", "\\1", exposure),
      exposure_study == "pietzner_2021_PMID34648354" ~ sub(".*_([^_]+_[^_]+)$", "\\1", exposure),
      TRUE ~ NA_character_
    ),
    UniProt = case_when(
      exposure_study == "sun_2023_PMID37794186" ~ sub("^[^_]+_([^_]+)_.*", "\\1", exposure),
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(SeqId = str_replace_all(SeqId, "_", "-")) %>%
  select(SeqId, UniProt, outcome, outcome_sex, method, BETA, SE, P, SNP, nsnp, 
         id.exposure, exposure_study, exposure_population, exposure_sex, exposure_cohort, exposure_finemap, exposure_samplesize, 
         id.outcome, outcome_study, outcome_population, outcome_sex, outcome_cohort, outcome_samplesize) %>%
  filter(
    method %in% c("Inverse variance weighted (multiplicative random effects)", "Wald ratio") &
      (method != "Wald ratio" | exposure_finemap == "coloc.abf"),
    exposure_finemap != "Finimom"
  ) %>%
  filter(!is.na(BETA)) %>%
  mutate(analysis_mr = "forward",
         ID = paste0(SeqId, UniProt, outcome, outcome_sex))

# data_mr reverse ====
LIST_FILES <- list.files("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/002_exposure-cancer_outcome-proteins/MR/", pattern = ".txt", all.files = T, full.names = T)
LIST_FILES <- LIST_FILES[str_detect(LIST_FILES, str_c(c(ID_seqid, ID_uniprot), collapse = "|"))]
data_mr <- lapply(LIST_FILES, data.table::fread, header = T, sep = "\t")
data_mr <- bind_rows(data_mr)
## format for seqid/uniprot matching
data_mr_reverse <- data_mr %>%
  mutate(
    SeqId = case_when(
      outcome_study == "ferkingstad_2021_PMID34857953" ~ sub("^([A-Za-z0-9]+_[A-Za-z0-9]+).*", "\\1", outcome),
      outcome_study == "zhang_2022_PMID35501419" ~ sub("^SeqId_(.*)", "\\1", outcome),
      outcome_study == "pietzner_2021_PMID34648354" ~ sub(".*_([^_]+_[^_]+)$", "\\1", outcome),
      TRUE ~ NA_character_
    ),
    UniProt = case_when(
      outcome_study == "sun_2023_PMID37794186" ~ sub("^[^_]+_([^_]+)_.*", "\\1", outcome),
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(SeqId = str_replace_all(SeqId, "_", "-")) %>%
  select(SeqId, UniProt, exposure, outcome_sex, method, BETA, SE, P, nsnp,
         id.exposure, exposure_study, exposure_population, exposure_sex, exposure_cohort, exposure_finemap, exposure_samplesize,
         id.outcome, outcome_study, outcome_population, outcome_sex, outcome_cohort, outcome_samplesize) %>%
  filter(method %in% c("Inverse variance weighted (multiplicative random effects)", "Wald ratio")) %>%
  filter(!is.na(BETA)) %>%
  filter(!str_detect(id.exposure, "NA")) %>%
  mutate(analysis_mr = "reverse",
         ID = paste0(SeqId, UniProt, exposure, exposure_sex)) %>%
  filter(ID %in% data_mr_forward$ID)

# format and filter for merging ====
df1 <- data_mr_forward %>%
  filter(exposure_finemap %in% c("SuSiE", "p_ld-0.001", "coloc.abf"),
         outcome_study == "huyghe_2018_PMID30510241",
         outcome_population == "ALL") %>%
  select(SeqId, UniProt, outcome, outcome_sex, analysis_mr, exposure_study, exposure_finemap, outcome_study, exposure_samplesize, outcome_samplesize, method, SNP, nsnp, BETA, SE, P) %>%
  rename(cancer = outcome,
         cancer_sex = outcome_sex)

df2 <- data_mr_reverse %>%
  filter(exposure_study == "huyghe_2018_PMID30510241",
         exposure_population == "ALL") %>%
  select(SeqId, UniProt, exposure, exposure_sex, analysis_mr, exposure_study, outcome_study, exposure_finemap, exposure_samplesize, outcome_samplesize, method, nsnp, BETA, SE, P) %>%
  rename(cancer = exposure,
         cancer_sex = exposure_sex)

data_mr <- bind_rows(df1, df2) 

## combine ====
table <- data %>%
  left_join(data_mr, by = c("SeqId", "cancer", "sex" = "cancer_sex")) %>%
  select(-UniProt.y) %>%
  rename(UniProt = UniProt.x) %>%
  droplevels() %>%
  mutate(
    cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
    sex = factor(sex, levels = c("combined", "female", "male"))
  ) %>%
  group_by(SeqId, Target, cancer, sex) %>%
  filter(all(c("forward", "reverse") %in% analysis_mr)) %>%
  ungroup()

df <- table %>%
  select(SeqId, UniProt, Target, TargetFullName, cancer, sex, coef, se_robust, FDR_BH) %>%
  rename(P = FDR_BH,
         BETA = coef,
         SE = se_robust) %>%
  mutate(exposure_study = "EPIC",
         method = "Cox") %>%
  unique()

df1 <- table %>%
  select(SeqId, UniProt, Target, TargetFullName, cancer, sex, BETA, SE, P, analysis_mr, SNP, nsnp, exposure_study, exposure_samplesize, outcome_study, outcome_samplesize, exposure_finemap, method) %>%
  unique()

table <- bind_rows(df,df1) %>%
  mutate(
    cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    exposure_finemap = factor(exposure_finemap, levels = c("coloc.abf", "p_ld-0.001", "SuSiE")),
    exposure_study = factor(exposure_study, levels = c("EPIC", "ferkingstad_2021_PMID34857953", "pietzner_2021_PMID34648354", "zhang_2022_PMID35501419", "huyghe_2018_PMID30510241"))
  ) %>%
  filter(sex == "combined") %>%
  arrange(Target) %>%
  select(SeqId, UniProt, Target, TargetFullName, cancer, sex, BETA, SE, P, analysis_mr, SNP, nsnp, method, exposure_study, exposure_finemap, exposure_samplesize, outcome_study, outcome_samplesize)

table <- table %>%
  filter(exposure_finemap != "SuSiE")

## save ====
write.table(x = table, file = "analysis/tables/manuscript/comparison-epic-MR.txt", 
            quote = FALSE, sep = "\t", row.names = F, col.names = T)

wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "comparison-epic-MR")
writeData(wb = wb, sheet = "comparison-epic-MR", table, colNames = T, rowNames = F)
saveWorkbook(wb, "analysis/tables/manuscript/comparison-epic-MR.xlsx", overwrite = TRUE)
