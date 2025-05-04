rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(openxlsx)

# cox results ====
data <- data.table::fread("analysis/002_coxph/003_format/001_analysis/cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP/formatted-analysis.txt")
data <- data %>%
  filter(cancer != "earlyonset",
         Organism == "Human") %>%
  mutate(cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
         sex = factor(sex, levels = c("combined", "female", "male")),
         followup = factor(followup, levels = c(0, 2, 5)),
         raw_complete = factor(raw_complete, levels = c("raw", "complete"))) %>%
  arrange(cancer, sex, followup, raw_complete, data_id, Target, SeqId, model_type)
## calculate FDR_BH within cancer, sex, raw_complete, followup, model_type, Organism; for Human proteins this will give an FDR_BH for 7,335 proteins
data <- data %>%
  group_by(cancer, sex, raw_complete, followup, model_type, Organism) %>%
  mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  ungroup() 
## data_plot
table(data$data_id)
data <- data %>%
  filter(
    raw_complete == "raw",
    model_type %in% c("model_1_extra", "model_2_extra")
  ) %>%
  select(SeqId, Target, TargetFullName, UniProt, cancer, sex, followup, model_type, n, nevent, 
         coef, coef_exp, se, se_robust, ci_lower_exp, ci_upper_exp, pval, FDR_BH, 
         degrees_freedom, concordance, concordance_se, likelihood_ratio_test, likelihood_ratio_test_pval, 
         wald_test, wald_test_pval, score_test, score_test_pval, robust_score_test, robust_score_test_pval)

# save ====
write.table(data, "analysis/tables/manuscript/results-epic.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "results-epic")
writeData(wb = wb, sheet = "results-epic", data, colNames = T, rowNames = F)
saveWorkbook(wb, "analysis/tables/manuscript/results-epic.xlsx", overwrite = TRUE)
