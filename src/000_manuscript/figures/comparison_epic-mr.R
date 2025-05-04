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

# MR comparison ====
## data
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
  select(SeqId, UniProt, Target, cancer, sex, coef, se_robust, pval, FDR_BH)
ID_seqid <- data$SeqId %>%
  unique() %>%
  str_replace_all("-", "_")
ID_uniprot <- unique(data$UniProt) 

## data_mr
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
data_mr <- data_mr %>%
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
  filter(!is.na(BETA))

## combine 
data_plot <- data %>%
  left_join(data_mr, by = "SeqId", suffix = c("", "_temp")) %>%
  mutate(SeqId = coalesce(SeqId, UniProt_temp)) %>% 
  left_join(data_mr, by = c("SeqId" = "UniProt"), suffix = c("", "_temp_2")) %>%
  select(-ends_with("_temp"), -ends_with("_temp_2")) %>%
  filter(cancer == outcome, sex == outcome_sex) %>%
  droplevels() %>%
  filter(exposure_finemap %in% c("SuSiE", "p_ld-0.001", "coloc.abf")) %>%
  filter(outcome_study == "huyghe_2018_PMID30510241") %>%
  filter(outcome_population == "ALL") %>%
  mutate(
    cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
    sex = factor(sex, levels = c("combined", "female", "male"))
  )

## forest plot
df <- data_plot %>%
  select(SeqId, Target, cancer, sex, coef, se_robust, FDR_BH) %>%
  rename(P = FDR_BH,
         BETA = coef,
         SE = se_robust) %>%
  mutate(exposure_study = "EPIC",
         exposure_finemap = "Cox")
df1 <- data_plot %>%
  select(SeqId, Target, cancer, sex, BETA, SE, P, exposure_study, exposure_finemap) 
df_plot <- bind_rows(df,df1) %>%
  mutate(
    cancer = factor(cancer, levels = c("overall", "colon", "rectum")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    exposure_finemap = factor(exposure_finemap, levels = c("Cox", "coloc.abf", "p_ld-0.001", "SuSiE")),
    exposure_study = factor(exposure_study, levels = c("EPIC", "ferkingstad_2021_PMID34857953", "pietzner_2021_PMID34648354", "zhang_2022_PMID35501419"))
  ) %>%
  filter(sex == "combined") %>%
  filter(exposure_finemap != "SuSiE") %>%
  unique() %>%
  arrange(Target)  # Sort Target alphabetically

## how many EPIC
dim(df_plot %>%
      filter(exposure_study == "EPIC") %>%
      select(Target, SeqId, cancer, sex) %>%
      unique())
dim(df_plot %>%
      filter(exposure_study == "EPIC") %>%
      select(Target, SeqId) %>%
      unique())
## how many SuSiE
dim(df_plot %>%
      filter(exposure_finemap == "SuSiE") %>%
      select(Target, SeqId, cancer, sex) %>%
      unique())
dim(df_plot %>%
      filter(exposure_finemap == "SuSiE") %>%
      select(Target, SeqId) %>%
      unique())
## how many coloc.abf
dim(df_plot %>%
      filter(exposure_finemap == "coloc.abf") %>%
      select(Target, SeqId, cancer, sex) %>%
      unique())
dim(df_plot %>%
      filter(exposure_finemap == "coloc.abf") %>%
      select(Target, SeqId) %>%
      unique())
## how many p_ld-0.001
dim(df_plot %>%
      filter(exposure_finemap == "p_ld-0.001") %>%
      select(Target, SeqId, cancer, sex) %>%
      unique())
dim(df_plot %>%
      filter(exposure_finemap == "p_ld-0.001") %>%
      select(Target, SeqId) %>%
      unique())

## plot all
list_plots <- lapply(split(df_plot, list(df_plot$cancer, df_plot$sex), drop = TRUE), function(sub_data) {
  cancer_level <- unique(sub_data$cancer)
  sex_level <- unique(sub_data$sex)
  
  p <- functions::forestplot(df = sub_data, 
                             name = Target, 
                             estimate = BETA, 
                             se = SE, 
                             pvalue = P, 
                             colour = exposure_study, 
                             shape = exposure_finemap, 
                             logodds = FALSE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    coord_cartesian(xlim = c(min(sub_data$BETA), max(sub_data$BETA))) +
    labs(title = paste(cancer_level),
         colour = "study", shape = "method") +
    theme(legend.position = "none")  # Remove legend from individual plots
  
  return(p)
})

# Extract legend from one of the plots (using the first non-null plot)
legend_plot <- functions::forestplot(df = df_plot, 
                                     name = Target, 
                                     estimate = BETA, 
                                     se = SE, 
                                     pvalue = P, 
                                     colour = exposure_study, 
                                     shape = exposure_finemap, 
                                     logodds = FALSE, 
                                     psignif = 0.05, 
                                     ci = 0.95) +
  labs(colour = "study", shape = "method") +
  theme(legend.position = "right")
legend <- cowplot::get_legend(legend_plot)  # Extract legend
# Arrange plots in a grid with the legend at the bottom
plot <- plot_grid(plotlist = list_plots, ncol = 2)
plot <- plot_grid(plot, legend, ncol = 2, rel_widths = c(1, 0.2))
# save
tiff(filename = "analysis/figures/manuscript/comparison-epic-mr.tiff",
     width = 1600, height = 800, units = "px")
plot
dev.off()
