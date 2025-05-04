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

# figure 1: volcano plot all ====
## data
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
    model_type == "model_2_extra",
    sex == "combined"
  )
## plot
p1 <- ggplot(data_plot, aes(x = coef_exp, y = -log10(FDR_BH))) +
  geom_point(aes(
    colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  scale_color_manual(values = c("YES" = palette[[2]][[3]], "NO" = "black")) +
  facet_grid(~cancer, 
             scales = "fixed") +
  theme_cowplot() +
  
  # lines 
  geom_vline(xintercept = 1, color = palette[[1]][[1]], linetype = "solid") +
  geom_hline(yintercept = -log10(0.05), color = palette[[1]][[1]], linetype = "solid") + 
  scale_alpha_continuous(range = c(0.1, 1)) +
  
  # labels
  ggrepel::geom_text_repel(
    data = data_plot %>%
      filter(FDR_BH < 0.05),  
    aes(x = coef_exp, 
        y = -log10(FDR_BH), 
        label = Target),  
    max.overlaps = 1000,
    min.segment.length = 0,
    colour = "black",
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5
  ) +
  
  # legend
  guides(alpha = "none", size = "none", colour = "none") +
  # labels
  xlab("Hazard ratio (95% confidence interval)") +
  ylab(label = "-log10(FDR pval)")
## save
tiff("analysis/figures/manuscript/results-epic.tiff", 
     width = 1600, height = 800, units = "px")
p1
dev.off()

