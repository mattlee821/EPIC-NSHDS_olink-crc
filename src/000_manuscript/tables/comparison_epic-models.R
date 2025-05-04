rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(openxlsx)

# model comparison ====
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
data_plot <- data %>%
  filter(
    followup == 0,
    raw_complete == "raw",
    model_type %in% c("model_1_extra", "model_2_extra")
  ) %>%
  select(exposure, cancer, sex, model_type, coef, FDR_BH)
## df combined
df1 <- data_plot %>%
  filter(model_type == "model_1_extra") %>%
  select(-model_type)
df2 <- data_plot %>%
  filter(model_type == "model_2_extra") %>%
  select(-model_type)
data_plot <- df1 %>%
  rename(coef_data1 = coef, 
         FDR_BH_data1 = FDR_BH) %>%
  left_join(
    df2 %>% rename(coef_data2 = coef, 
                   FDR_BH_data2 = FDR_BH),
    by = c("exposure", "cancer", "sex")
  ) %>%
  mutate(
    association = case_when(
      FDR_BH_data1 < 0.05 & FDR_BH_data2 < 0.05 ~ "both",   # Both significant
      FDR_BH_data1 < 0.05 ~ "model 1",                     # Only FDR_BH_data1 significant
      FDR_BH_data2 < 0.05 ~ "model 2",                       # Only FDR_BH_data2 significant
      TRUE ~ "none"                                      # All other points
    )
  )
## cor
table <- data_plot %>%
  group_by(cancer, sex) %>%
  summarise(
    cor_all = cor(coef_data1, coef_data2, use = "complete.obs", method = "pearson"),  # Full data correlation
    cor_fdr = cor(coef_data1[FDR_BH_data2 < 0.05], coef_data2[FDR_BH_data2 < 0.05], use = "complete.obs", method = "pearson"),  # Significant data correlation
    .groups = "drop"
  ) %>%
  select(cancer, sex, cor_all, cor_fdr)
## save
write.table(x = table, file = "analysis/tables/manuscript/comparison-epic-models.txt", 
            quote = FALSE, sep = "\t", row.names = F, col.names = T)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "comparison-epic-models")
writeData(wb = wb, sheet = "comparison-epic-models", table, colNames = T, rowNames = F)
saveWorkbook(wb, "analysis/tables/manuscript/comparison-epic-models.xlsx", overwrite = TRUE)
