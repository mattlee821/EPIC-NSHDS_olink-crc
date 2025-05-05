rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)

# MR comparison ====
## data
main_file <- "analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-FALSE_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt"
temp <- data.table::fread("data/processed/EPIC_olink-immuneonc_feature-metadata.txt") %>%
  dplyr::select(OlinkID, UNIPROT, Assay)
data <- data.table::fread(main_file)
data <- data %>%
  mutate(
    outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    followup = factor(followup, levels = c(0, 2, 5))
  ) %>%
  dplyr::arrange(exposure, outcome, sex, followup, model) %>%
  dplyr::group_by(exposure, outcome, sex, followup, model) %>%
  dplyr::mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(followup == 0, model == "model_2") %>%
  dplyr::left_join(temp, by = c("exposure" = "OlinkID")) %>%
  dplyr::rename(UniProt = UNIPROT)

ID_seqid <- NULL
ID_uniprot <- unique(data$UniProt) 

## data_mr
LIST_FILES <- list.files("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/001_exposure-proteins_outcome-cancer/MR/", pattern = ".txt", all.files = T, full.names = T)
LIST_FILES <- LIST_FILES[stringr::str_detect(LIST_FILES, stringr::str_c(c(ID_seqid, ID_uniprot), collapse = "|"))]
data_mr <- lapply(LIST_FILES, data.table::fread, header = T, sep = "\t")
data_mr <- dplyr::bind_rows(data_mr)
## remove all coloc.abf 
data_mr <- data_mr %>%
  dplyr::filter(!(exposure_finemap == "coloc.abf")) 
## add in all singlesnp MR for coloc.abf
LIST_FILES <- list.files("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/001_exposure-proteins_outcome-cancer/singlesnp/", pattern = ".txt", all.files = T, full.names = T)
LIST_FILES <- LIST_FILES[stringr::str_detect(LIST_FILES, stringr::str_c(c(ID_seqid, ID_uniprot), collapse = "|"))]
data_mrsinglesnp <- lapply(LIST_FILES, data.table::fread, header = T, sep = "\t")
data_mrsinglesnp <- dplyr::bind_rows(data_mrsinglesnp)
data_mrsinglesnp <- data_mrsinglesnp %>%
  dplyr::filter(exposure_finemap == "coloc.abf") %>%
  dplyr::mutate(method = "Wald ratio")
data_mr <- data_mr %>%
  dplyr::bind_rows(data_mrsinglesnp)
## format for seqid/uniprot matching
data_mr <- data_mr %>%
  dplyr::mutate(
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
  dplyr::mutate(SeqId = stringr::str_replace_all(SeqId, "_", "-")) %>%
  dplyr::select(SeqId, UniProt, outcome, outcome_sex, method, BETA, SE, P, SNP, nsnp, 
         id.exposure, exposure_study, exposure_population, exposure_sex, exposure_cohort, exposure_finemap, exposure_samplesize, 
         id.outcome, outcome_study, outcome_population, outcome_sex, outcome_cohort, outcome_samplesize) %>%
  dplyr::filter(
    method %in% c("Inverse variance weighted (multiplicative random effects)", "Wald ratio") &
      (method != "Wald ratio" | exposure_finemap == "coloc.abf"),
    exposure_finemap != "Finimom"
  ) %>%
  dplyr::filter(!is.na(BETA))

## combine 
data_plot <- data %>%
  dplyr::left_join(data_mr, by = c("UniProt", "outcome", "sex" = "outcome_sex")) %>%
  dplyr::filter(exposure_finemap %in% c("p_ld-0.001", "coloc.abf")) %>%
  dplyr::filter(outcome_study == "huyghe_2018_PMID30510241") %>%
  dplyr::filter(outcome_population == "ALL") %>%
  dplyr::filter(exposure_population == "ALL") %>%
  dplyr::mutate(
    cancer = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male"))
  )

temp <- data_plot %>%
  dplyr::filter(Assay == "IL8",
                cancer == "overall", 
                sex == "combined",
                exposure_population == "ALL",
                outcome_population == "ALL")

## forest plot
df <- data_plot %>%
  select(Assay, UniProt, cancer, sex, coef, se, FDR_BH) %>%
  rename(P = FDR_BH,
         BETA = coef,
         SE = se) %>%
  mutate(exposure_study = "EPIC",
         exposure_finemap = "log")
df1 <- data_plot %>%
  select(Assay, UniProt, cancer, sex, BETA, SE, P, exposure_study, exposure_finemap) 
df_plot <- bind_rows(df,df1) %>%
  mutate(
    cancer = factor(cancer, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    exposure_finemap = factor(exposure_finemap, levels = c("log", "coloc.abf", "p_ld-0.001", "SuSiE"))
    # exposure_study = factor(exposure_study, levels = c("EPIC", "ferkingstad_2021_PMID34857953", "pietzner_2021_PMID34648354", "zhang_2022_PMID35501419"))
  ) %>%
  filter(sex == "combined") %>%
  filter(exposure_finemap != "SuSiE") %>%
  unique() %>%
  group_by(UniProt, exposure_study, exposure_finemap, cancer, sex) %>%
  mutate(row_index = case_when(
    exposure_study == "sun_2023_PMID37794186" ~ paste0("OlinkID: ", as.character(row_number())),
    exposure_study == "EPIC" ~ "EPIC",
    TRUE ~ NA_character_
  )) %>%
  ungroup() %>%
  arrange(Assay)  # Sort Assay alphabetically

## plot all
list_plots <- lapply(split(df_plot, list(df_plot$cancer, df_plot$sex), drop = TRUE), function(sub_data) {
  cancer_level <- unique(sub_data$cancer)
  sex_level <- unique(sub_data$sex)
  
  p <- functions::forestplot(df = sub_data, 
                             name = Assay, 
                             estimate = BETA, 
                             se = SE, 
                             pvalue = P, 
                             colour = row_index, 
                             shape = exposure_finemap, 
                             logodds = FALSE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    ggplot2::coord_cartesian(xlim = c(min(sub_data$BETA), max(sub_data$BETA))) +
    ggplot2::labs(title = paste(cancer_level),
         colour = "study", shape = "method") +
    ggplot2::theme(legend.position = "none")  # Remove legend from individual plots
  
  return(p)
})

# Extract legend from one of the plots (using the first non-null plot)
legend_plot <- functions::forestplot(df = df_plot, 
                                     name = Assay, 
                                     estimate = BETA, 
                                     se = SE, 
                                     pvalue = P, 
                                     colour = row_index, 
                                     shape = exposure_finemap, 
                                     logodds = FALSE, 
                                     psignif = 0.05, 
                                     ci = 0.95) +
  ggplot2::labs(colour = "study", shape = "method") +
  ggplot2::theme(legend.position = "right")
legend <- cowplot::get_legend(legend_plot)  # Extract legend
# Arrange plots in a grid with the legend at the bottom
plot <- cowplot::plot_grid(plotlist = list_plots, ncol = 3)
plot <- cowplot::plot_grid(plot, legend, ncol = 2, rel_widths = c(1, 0.1))
# save
tiff(filename = "manuscript/figures/comparison_epic-immuneonc-mr.tiff",
     width = 1600, height = 1000, units = "px")
plot
dev.off()
