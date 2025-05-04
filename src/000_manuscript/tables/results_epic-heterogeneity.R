rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(openxlsx)

# data_epic ====
data <- data.table::fread("analysis/002_coxph/003_format/001_analysis/cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP/formatted-analysis.txt")
data_epic <- data %>%
  filter(cancer != "earlyonset",
         Organism == "Human") %>%
  mutate(cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
         sex = factor(sex, levels = c("combined", "female", "male")),
         followup = factor(followup, levels = c(0, 2, 5)),
         raw_complete = factor(raw_complete, levels = c("raw", "complete"))) %>%
  arrange(cancer, sex, followup, raw_complete, data_id, Target, SeqId, model_type)
## calculate FDR_BH within cancer, sex, raw_complete, followup, model_type, Organism; for Human proteins this will give an FDR_BH for 7,335 proteins
data_epic <- data_epic %>%
  group_by(cancer, sex, raw_complete, followup, model_type, Organism) %>%
  mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  ungroup()
## filter
data_epic <- data_epic %>%
  filter(
    followup == 0,
    raw_complete == "raw",
    model_type == "model_2_extra",
    FDR_BH < 0.05)

id_sex <- data_epic %>%
  mutate(id = paste0(exposure, cancer, sex)) %>%
  filter(sex == "combined") %>%
  pull(id) %>%
  unique()

id_site <- data_epic %>%
  mutate(id = paste0(exposure, cancer, sex)) %>%
  filter(cancer == "overall") %>%
  pull(id) %>%
  unique()

## data heterogeneity sex ====
data_heterogeneity <- data.table::fread("analysis/002_coxph/003_format/002_heterogeneity/sex/cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP/formatted-heterogeneity.txt")
data_heterogeneity <- data_heterogeneity %>%
  filter(cancer != "earlyonset",
         Organism == "Human") %>%
  mutate(cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
         sex = factor(sex, levels = c("combined", "female", "male")),
         followup = factor(followup, levels = c(0, 2, 5)),
         raw_complete = factor(raw_complete, levels = c("raw", "complete"))) %>%
  arrange(cancer, sex, followup, raw_complete, data_id, Target, SeqId, group)
## filter
data_heterogeneity_sex <- data_heterogeneity %>%
  filter(
    followup == 0,
    sex == "combined",
    raw_complete == "raw",
    group == "model2_extra") %>%
  mutate(id = paste0(exposure, cancer, sex)) %>%
  filter(id %in% id_sex) %>%
  group_by(cancer, sex, raw_complete, followup, group, Organism) %>%
  mutate(FDR_BH = p.adjust(p = pval_heterogeneity, method = "BH")) %>%
  ungroup() %>%
  select(SeqId, UniProt, Target, TargetFullName, cancer, sex, pval_heterogeneity)

## data heterogeneity site ====
data_heterogeneity <- data.table::fread("analysis/002_coxph/003_format/002_heterogeneity/cancer/cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP/formatted-heterogeneity.txt")
data_heterogeneity <- data_heterogeneity %>%
  filter(Organism == "Human") %>%
  mutate(sex = factor(sex, levels = c("combined", "female", "male")),
         followup = factor(followup, levels = c(0, 2, 5)),
         raw_complete = factor(raw_complete, levels = c("raw", "complete"))) %>%
  arrange(cancer, sex, followup, raw_complete, data_id, Target, model_type, SeqId)
## filter
data_heterogeneity_site <- data_heterogeneity %>%
  mutate(id = paste0(exposure, cancer, sex)) %>%
  filter(
    followup == 0,
    raw_complete == "raw",
    model_type == "model_2_heterogeneity_extra",
    id %in% id_site
  ) %>%
  group_by(cancer, sex, raw_complete, followup, model_type, Organism) %>%
  mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  ungroup() %>%
  select(SeqId, UniProt, Target, TargetFullName, cancer, sex, pval)

# table ====
table <- full_join(
  data_heterogeneity_sex %>%
    select(SeqId, UniProt, Target, TargetFullName, cancer, sex, pval_heterogeneity),  
  data_heterogeneity_site %>%
    select(SeqId, UniProt, Target, TargetFullName, cancer, sex, pval),  
  by = c("SeqId", "UniProt", "Target", "TargetFullName", "cancer", "sex")
) %>%
  rename(
    pval_heterogeneity_sex = pval_heterogeneity,
    pval_heterogeneity_site = pval) %>%
  select(SeqId, UniProt, Target, TargetFullName, cancer, sex, pval_heterogeneity_sex, pval_heterogeneity_site)

## save ====
write.table(x = table, file = "analysis/tables/manuscript/results-epic-heterogeneity.txt", 
            quote = FALSE, sep = "\t", row.names = F, col.names = T)

wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "epic-heterogeneity")
writeData(wb = wb, sheet = "epic-heterogeneity", table, colNames = T, rowNames = F)
saveWorkbook(wb, "analysis/tables/manuscript/results-epic-heterogeneity.xlsx", overwrite = TRUE)
