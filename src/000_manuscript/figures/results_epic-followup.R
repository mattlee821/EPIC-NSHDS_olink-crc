rm(list=ls())
set.seed(821)

# environment ====

# data ====
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
id_associations <- data %>%
  filter(
    followup == 0,
    raw_complete == "raw",
    model_type == "model_2_extra",
    FDR_BH < 0.05
  ) %>%
  select(exposure, cancer, sex) %>%
  mutate(id = paste(exposure, ";", cancer, ";", sex)) %>%
  pull(id)

## data_plot
data_plot <- data %>%
  filter(
    raw_complete == "raw",
    model_type == "model_2_extra",
  ) %>%
  mutate(id = paste(exposure, ";", cancer, ";", sex)) %>%
  filter(id %in% id_associations,
         sex == "combined")

## plot ====
list_plots <- lapply(split(data_plot, list(data_plot$cancer, data_plot$sex), drop = TRUE), function(sub_data) {
  cancer_level <- unique(sub_data$cancer)
  sex_level <- unique(sub_data$sex)
  
  p <- functions::forestplot(df = sub_data, 
                             name = Target, 
                             estimate = coef, 
                             se = se_robust, 
                             pvalue = FDR_BH, 
                             colour = followup, 
                             logodds = T, 
                             psignif = 0.05, 
                             ci = 0.95) +
    labs(title = paste(cancer_level),
         colour = "study", shape = "method") +
    xlab("Hazard ratio (95% confidence interval)") +
    theme(legend.position = "none")  # Remove legend from individual plots
  
  return(p)
})

# Extract legend from one of the plots (using the first non-null plot)
legend_plot <- functions::forestplot(df = data_plot, 
                                     name = Target, 
                                     estimate = coef, 
                                     se = se_robust, 
                                     pvalue = FDR_BH, 
                                     colour = followup, 
                                     logodds = T, 
                                     psignif = 0.05, 
                                     ci = 0.95) +
  labs(colour = "followup exclusion") +
  xlab("Hazard ratio (95% confidence interval)") +
  theme(legend.position = "right")
legend <- cowplot::get_legend(legend_plot)  # Extract legend
# Arrange plots in a grid with the legend at the bottom
plot <- plot_grid(plotlist = list_plots, ncol = 3)
plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(1, 0.2))

## save
tiff("analysis/figures/manuscript/results-epic-followup.tiff", 
     width = 1600, height = 800, units = "px")
plot
dev.off()
