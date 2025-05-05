rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(dplyr)
library(SomaDataIO)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggrepel)

# protein info ====
info_analyte <- read_adat("/data/Epic/subprojects/Depot_Somalogic/sources/SS-2230719.SS-2215915.SS-2230718.SS-2230720_v4.1_CitratePlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
info_analyte <- getAnalyteInfo(info_analyte)
info_analyte <- info_analyte %>%
  select(SeqId, UniProt, Target)

# data ====
LIST_FILES <- list.files(path = "analysis/002_coxph/002_formatted/associations/", pattern = ".txt", all.files = T, full.names = T)
data <- lapply(X = LIST_FILES, fread, colClasses = c("EntrezGeneID" = "character"))
data <- bind_rows(data)
data <- data %>%
  filter(followup == 0,
         model_type == 2)
data <- data %>%
  select(SeqId, UniProt, Target, cancer, sex, coef, se, pval) %>%
  unique()
id_UniProt <- unique(data$UniProt)
id_SeqId <- unique(data$SeqId)

# data_MR ====
LIST_FILES <- c(
  list.files(path = "/data/MET_share/work/001_projects/adiposity-proteins-crc/analysis/002_colocalization/004_coloc/125k", pattern = ".txt", all.files = T, full.names = T),
  list.files(path = "/data/MET_share/work/001_projects/adiposity-proteins-crc/analysis/002_colocalization/004_coloc/250k", pattern = ".txt", all.files = T, full.names = T),
  list.files(path = "/data/MET_share/work/001_projects/adiposity-proteins-crc/analysis/002_colocalization/004_coloc/500k", pattern = ".txt", all.files = T, full.names = T),
  list.files(path = "/data/MET_share/work/001_projects/adiposity-proteins-crc/analysis/002_colocalization/004_coloc/1mb/", pattern = ".txt", all.files = T, full.names = T)
)
data_coloc <- lapply(X = LIST_FILES, fread)
data_coloc <- bind_rows(data_coloc)

# format ====
data_coloc <- data_coloc %>%
  select(-exposure) %>%
  separate(id, into = c("exposure", "SNP", "window", "study_exposure"), sep = ";") %>%
  mutate(outcome = gsub(pattern = ".txt", replacement = "", x = outcome)) %>%
  separate(outcome, into = c("study_outcome", "outcome"), sep = ";") %>%
  separate(outcome, into = c("cancer", "sex", "cohort_outcome", "ancestry"), sep = "_") %>%
  mutate(
    ancestry = if_else(is.na(ancestry), "ALL", ancestry),
    ancestry = if_else(ancestry == "allancestries", "ALL", ancestry),
    ancestry = if_else(ancestry == "european", "EU", ancestry),
    cancer = if_else(cancer == "crc", "overall", cancer), 
    cancer = if_else(cancer == "rectal", "rectum", cancer),
    cohort_outcome = ifelse(is.na(cohort_outcome), "CRC_early-onset", cohort_outcome)
  )

# filter ====
data_coloc <- data_coloc %>%
  filter(study_outcome == "huyghe_2018_PMID30510241") %>%
  mutate(window = case_when(
    window == "1mb" ~ "2mb",
    window == "500k" ~ "1mb",
    window == "250k" ~ "500k",
    window == "125k" ~ "250k",
    TRUE ~ window
  )) %>%
  mutate(
    SeqId = if_else(population == "DECODE", str_replace(str_extract(exposure, "^[^_]*_[^_]*"), "_", "-"), NA_character_),
    UniProt = if_else(population != "DECODE", str_extract(exposure, "(?<=_)[^_]*(?=_[^_]*$)"), NA_character_)
  )

# combine with protein info ====
a <- data_coloc %>%
  filter(population == "DECODE") %>%
  left_join(info_analyte, by = "SeqId") %>%
  select(-UniProt.x) %>%
  rename(UniProt = UniProt.y)
b <- data_coloc %>%
  filter(population != "DECODE") %>%
  left_join(info_analyte, by = "UniProt", relationship = "many-to-many") %>%
  select(-SeqId.x) %>%
  rename(SeqId = SeqId.y)
data_coloc <- bind_rows(a,b)

# filter for associations ====
data_coloc <- data_coloc %>%
  filter(SeqId %in% id_SeqId | UniProt %in% id_UniProt)
length(unique(data_coloc$SeqId))
length(unique(data_coloc$UniProt))
length(unique(data_coloc$Target))

# filter for combinations ====
data$id <- paste0(data$SeqId, ";", data$cancer, ";", data$sex)
data_coloc$id <- paste0(data_coloc$SeqId, ";", data_coloc$cancer, ";", data_coloc$sex)
data_coloc <- data_coloc %>%
  filter(id %in% data$id) %>%
  filter(ancestry != "EU") %>%
  select(SeqId, UniProt, Target, cancer, ancestry, population, sex, SNP, window, nsnps, h0,h1,h2,h3,h4, prior_p1,prior_p2,prior_p12)
length(unique(data_coloc$SeqId))
length(unique(data_coloc$UniProt))
length(unique(data_coloc$Target))
data_coloc$priors <- paste0(data_coloc$prior_p1, ";", data_coloc$prior_p2, ";", data_coloc$prior_p12)

# associations ====
data_associations <- data_coloc %>%
  filter(h4 > 0.8 & window == "1mb" & prior_p1 == 1E-6 & prior_p2 == 1E-6 & prior_p12 == 1E-7)

# table ====
table <- data_coloc %>%
  select(SeqId, UniProt, Target, population, cancer, sex, SNP, window, nsnps, h0,h1,h2,h3,h4, priors) %>%
  rename(study_protein = population) %>%
  arrange(Target, cancer, sex) 
write.table(x = table, file = "analysis/004_colocalization/table_associations.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# plot ====
plot_data <- data_coloc %>%
  mutate(
    cancer = factor(cancer, levels = c("overall", "colon", "rectum", "earlyonset")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    window = factor(window, levels = c("250k", "500k", "1mb", "2mb")),
    priors = factor(priors, levels = c("1e-04;1e-04;1e-05", "1e-04;1e-06;1e-07", "1e-06;1e-04;1e-07", "1e-06;1e-06;1e-07")),
    label = paste0(Target, "_", SeqId, "_", UniProt)
  ) %>%
  arrange(Target) 

p1 <- ggplot(data = plot_data, 
             aes(x = h0, y = h4, 
                 colour = window, shape = priors)) +
  geom_point() +
  theme_cowplot() +
  theme(legend.position = "none") +
  geom_label_repel(data = subset(plot_data, h4 > 0.8), 
                   aes(label = paste(Target, cancer, sex)), color = "black", min.segment.length = 0)

p2 <- ggplot(data = plot_data, 
             aes(x = h1, y = h4,
                 colour = window, shape = priors)) +
  geom_point() +
  theme_cowplot() +
  theme(legend.position = "none") +
  geom_label_repel(data = subset(plot_data, h4 > 0.8), 
                   aes(label = paste(Target, cancer, sex)), color = "black", min.segment.length = 0)

p3 <- ggplot(data = plot_data, 
             aes(x = h2, y = h4,
                 colour = window, shape = priors)) +
  geom_point() +
  theme_cowplot() +
  theme(legend.position = "none") +
  geom_label_repel(data = subset(plot_data, h4 > 0.8), 
                   aes(label = paste(Target, cancer, sex)), color = "black", min.segment.length = 0)

p4 <- ggplot(data = plot_data, 
             aes(x = h3, y = h4,
                 colour = window, shape = priors)) +
  geom_point() +
  theme_cowplot() +
  theme(legend.position = "none") +
  geom_label_repel(data = subset(plot_data, h4 > 0.8), 
                   aes(label = paste(Target, cancer, sex)), color = "black", min.segment.length = 0)

legend_plot <- ggplot(data = plot_data, 
                      aes(x = h3, y = h4,
                          colour = window, shape = priors)) +
  geom_point() +
  theme_cowplot() 
legend_plot <- cowplot::get_legend(legend_plot)
p <- plot_grid(... = p1,p2,p3,p4,legend_plot, nrow = 1, rel_widths = c(1,1,1,1,.4))

tiff(filename = "analysis/figures/plot_colocalisation_associations/colocalisation.tiff",
     height = 600, width = 1400, units = "px")    
p
dev.off()




