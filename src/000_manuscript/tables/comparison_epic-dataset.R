rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(openxlsx)
library(ggplot2)

# source ====
palette <- palette()

# dataset comparison ====
## data
list_files <- list.files(path = "analysis/002_coxph/003_format/001_analysis/", pattern = "formatted-analysis.txt", all.files = T, recursive = T, full.names = T)
main_file <- "analysis/002_coxph/003_format/001_analysis//cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP/formatted-analysis.txt"
data <- data.table::fread(main_file)
data <- data %>%
  filter(cancer != "earlyonset", Organism == "Human") %>%
  mutate(
    cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    followup = factor(followup, levels = c(0, 2, 5)),
    raw_complete = factor(raw_complete, levels = c("raw", "complete"))
  ) %>%
  arrange(cancer, sex, followup, raw_complete, data_id, Target, SeqId, model_type) %>%
  group_by(cancer, sex, raw_complete, followup, model_type, Organism) %>%
  mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  ungroup() %>%
  filter(followup == 0, raw_complete == "raw", model_type == "model_2_extra")
## table/plot loop
table <- data.frame(data = character(), cancer = character(), sex = character(), cor_all = numeric(), cor_fdr = numeric())
extract_second_last_component <- function(filepath) {
  path_components <- strsplit(filepath, "/")[[1]]
  second_last <- path_components[length(path_components) - 1]
  return(second_last)
}
# Process each file sequentially ====
for (file in list_files) {
  # Skip the main file
  if (file == main_file) next
  label <- extract_second_last_component(file)
  cat("Data:", label, "\n")
  # Read and format the new data
  df <- data.table::fread(file) 
  df <- df %>%
    filter(cancer != "earlyonset", Organism == "Human") %>%
    mutate(
      cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
      sex = factor(sex, levels = c("combined", "female", "male")),
      followup = factor(followup, levels = c(0, 2, 5)),
      raw_complete = factor(raw_complete, levels = c("raw", "complete"))
    ) %>%
    arrange(cancer, sex, followup, raw_complete, data_id, Target, SeqId, model_type) %>%
    group_by(cancer, sex, raw_complete, followup, model_type, Organism) %>%
    mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
    ungroup() %>%
    filter(followup == 0, raw_complete == "raw", model_type == "model_2_extra")
  
  # combined
  data_plot <- data %>%
    select(exposure, cancer, sex, coef_data1 = coef, FDR_BH_data1 = FDR_BH) %>%
    inner_join(
      df %>% select(exposure, cancer, sex, coef_data2 = coef, FDR_BH_data2 = FDR_BH),
      by = c("exposure", "cancer", "sex")
    ) %>%
    mutate(
      association = case_when(
        FDR_BH_data1 < 0.05 & FDR_BH_data2 < 0.05 ~ "both",   # Both significant
        FDR_BH_data1 < 0.05 ~ "main",                     # Only FDR_BH_data1 significant
        FDR_BH_data2 < 0.05 ~ "sensitivity",                       # Only FDR_BH_data2 significant
        TRUE ~ "none"                                      # All other points
      )
    )
  
  # Compute correlation between coef columns
  table_cor <- data_plot %>%
    group_by(cancer, sex) %>%
    summarise(
      cor_all = cor(coef_data1, coef_data2, use = "complete.obs", method = "pearson"),  # Full data correlation
      cor_fdr = cor(coef_data1[FDR_BH_data1 < 0.05], coef_data2[FDR_BH_data1 < 0.05], use = "complete.obs", method = "pearson"),  # Significant data correlation
      .groups = "drop"
    ) %>%
    mutate(data = extract_second_last_component(file)) %>%
    select(data, cancer, sex, cor_all, cor_fdr)
  table <- bind_rows(table, table_cor)
  cat("cor done \n")
  
  # plot
  p1 <- ggplot(data_plot, aes(x = coef_data2, y = coef_data1, colour = association)) +
    geom_point(aes(alpha = ifelse(association == "none", 0.5, 1))) +  
    scale_color_manual(
      values = c(
        "both" = "green",   
        "main" = "orange",  
        "sensitivity" = "blue",
        "none" = "black"    
      )) +
    geom_vline(xintercept = 0, color = palette[[1]][[1]], linetype = "solid") +
    geom_hline(yintercept = 0, color = palette[[1]][[1]], linetype = "solid") + 
    facet_grid(sex~cancer, 
               scales = "fixed") +
    labs(
      x = label,
      y = "Main analysis data"
    ) +
    cowplot::theme_cowplot() +
    guides(alpha = "none")  
  tiff(paste0("analysis/figures/manuscript/epic-data-comparison/correlation_", label, ".tiff"), 
       width = 20, height = 20, units = "cm", res = 100)
  print(p1)
  dev.off()
  cat("plot done \n")
  
}
write.table(x = table, file = "analysis/tables/manuscript/comparison-epic-data.txt", 
            quote = FALSE, sep = "\t", row.names = F, col.names = T)
wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "comparison-epic-data")
writeData(wb = wb, sheet = "comparison-epic-data", table, colNames = T, rowNames = F)
saveWorkbook(wb, "analysis/tables/manuscript/comparison-epic-data.xlsx", overwrite = TRUE)
