rm(list=ls())
set.seed(821)

# environment ====
library(functions)
library(ggplot2)
library(wesanderson)
library(dplyr)
library(cowplot)
library(data.table)
library(ggrepel)
library(tidyr)
library(data.table)
library(dplyr)
library(tidyr)
library(pathfindR)
library(cowplot)
library(grid)
library(UpSetR)
library(ggplot2)

# source ====
palette <- palette()

# heterogeneity by sex/site for overall/combined ====
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
## filter ====
data_plot <- data %>%
  filter(
    followup == 0,
    raw_complete == "raw",
    model_type == "model_2_extra"
  ) %>%
  select(exposure, Target, cancer, sex, coef_exp, se_robust, FDR_BH)
VARS <- data_plot %>%
  filter(FDR_BH < 0.05,
         sex == "combined") %>%
  pull(exposure)
## data heterogeneity ====
data_epic <- data.table::fread("analysis/tables/manuscript/supplement-table2_results-epic-cox.txt")
data_epic <- data_epic %>%
  filter(FDR_BH < 0.05,
         followup == 0,
         model_type == "model_2_extra",
         cancer == "overall",
         sex == "combined")
data_sex <- data.table::fread("analysis/tables/manuscript/supplement-table3_results-epic-heterogeneity.txt")
data_sex <- data_sex %>%
  filter(
    pval_heterogeneity_sex < 0.05,
    cancer == "overall",
    Target %in% data_epic$Target)
data_site <- data.table::fread("analysis/tables/manuscript/supplement-table3_results-epic-heterogeneity.txt")
data_site <- data_site %>%
  filter(
    pval_heterogeneity_site < 0.05,
    sex == "combined",
    Target %in% data_epic$Target)
id_sex <- unique(data_sex$Target)
id_site <- unique(data_site$Target)
id <- intersect(id_sex, id_site)

## format data_plot for labels ====
df1 <- data_plot %>%
  filter(Target %in% id_sex,
         !Target %in% id) %>%
  mutate(label = Target,
         heterogeneity = "sex")
df2 <- data_plot %>%
  filter(Target %in% id_site,
         !Target %in% id) %>%
  mutate(label = Target,
         heterogeneity = "site")
df3 <- data_plot %>%
  filter(Target %in% id) %>%
  mutate(label = Target,
         heterogeneity = "both")
df4 <- data_plot %>%
  filter(!Target %in% c(id_sex, id_site))
data_plot <- bind_rows(df1,df2,df3,df4)

## plot
p1 <- ggplot() +
  
  # First set of points where heterogeneity is NA
  geom_point(
    data = data_plot %>%
      filter(is.na(heterogeneity)), 
    aes(x = coef_exp, 
        y = -log10(FDR_BH),
        colour = ifelse(FDR_BH < 0.05, "YES", "NO"),
        alpha = ifelse(FDR_BH < 0.05, 1, 0.5))
  ) +
  
  # Second set of points where heterogeneity is NOT NA
  geom_point(
    data = data_plot %>%
      filter(!is.na(heterogeneity)),  
    aes(x = coef_exp, 
        y = -log10(FDR_BH),
        size = 2,
        fill = heterogeneity,
        colour = heterogeneity
    ),
    shape = 21
  ) +
  
  # Lines
  geom_vline(xintercept = 1, color = palette[[1]][[1]], linetype = "solid") +
  geom_hline(yintercept = -log10(0.05), color = palette[[1]][[1]], linetype = "solid") + 
  scale_alpha_continuous(range = c(0.1, 1)) +
  
  # Labels with correct dataset
  ggrepel::geom_text_repel(
    data = data_plot,
    aes(x = coef_exp, 
        y = -log10(FDR_BH), 
        label = label),
    max.overlaps = 1000,
    min.segment.length = 0, force = 10,
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5,
    colour = "black"
  ) +
  
  # Legend
  scale_color_manual(values = c("YES" = palette[[2]][[3]], "NO" = "black")) +
  guides(alpha = "none", colour = "none", size = "none") +
  xlab("Hazard ratio (95% confidence interval)") +
  ylab(label = "-log10(FDR pval)") +
  
  # Theme
  facet_grid(sex ~ cancer, scales = "fixed") +
  theme_cowplot() 

## save
tiff("analysis/figures/manuscript/results-epic-heterogeneity.tiff", 
     width = 1600, height = 800, units = "px")
p1
dev.off()
