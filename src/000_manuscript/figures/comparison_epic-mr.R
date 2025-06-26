rm(list=ls())
set.seed(821)

# figure: forestplot compare EPIC-MR

# environment ====
source("src/000_source.R")
rm(list=ls())

# data epic ====
data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt")
data_info <- data.table::fread("data/processed/EPIC_olink-explore.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% data$exposure)
data_explore <- data %>%
  dplyr::left_join(data_info, by = c("exposure"="OlinkID")) 

data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt")
data_info <- data.table::fread("data/processed/EPIC_olink-immuneonc.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% data$exposure)
data_immuneonc <- data %>%
  dplyr::left_join(data_info, by = c("exposure"="OlinkID"))

ID <- unique(c(data_explore$UniProt, data_immuneonc$UniProt))

# data_mr ====
LIST_FILES <- list.files("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/001_exposure-proteins_outcome-cancer/MR/", pattern = ".txt", all.files = T, full.names = T)
LIST_FILES <- LIST_FILES[stringr::str_detect(LIST_FILES, stringr::str_c(c(ID), collapse = "|"))]
data_mr <- lapply(LIST_FILES, data.table::fread, header = T, sep = "\t")
data_mr <- dplyr::bind_rows(data_mr)
## remove all coloc.abf 
data_mr <- data_mr %>%
  dplyr::filter(!(exposure_finemap == "coloc.abf")) 
## add in all singlesnp MR for coloc.abf
LIST_FILES <- list.files("/data/MET_share/work/001_projects/proteins-crc/analysis/001_MR/004_MR/001_exposure-proteins_outcome-cancer/singlesnp/", pattern = ".txt", all.files = T, full.names = T)
LIST_FILES <- LIST_FILES[stringr::str_detect(LIST_FILES, stringr::str_c(c(ID), collapse = "|"))]
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

# data_plot_immuneonc ====
data_plot <- data_immuneonc %>%
  dplyr::mutate(
    outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    followup = factor(followup, levels = c(0, 2, 5)),
    coef_exp = as.numeric(coef_exp)
  ) %>%
  dplyr::arrange(exposure, outcome, sex, followup, model) %>%
  dplyr::group_by(exposure, outcome, sex, followup, model) %>%
  dplyr::mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(followup == 0, model == "model_2") %>%
  dplyr::left_join(data_mr, by = c("UniProt", "outcome", "sex" = "outcome_sex")) %>%
  dplyr::filter(exposure_finemap %in% c("p_ld-0.001", "coloc.abf")) %>%
  dplyr::filter(outcome_study == "huyghe_2018_PMID30510241") %>%
  dplyr::filter(outcome_population == "ALL") %>%
  dplyr::filter(exposure_population == "ALL") %>%
  dplyr::mutate(
    cancer = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male"))
  )

df <- data_plot %>%
  dplyr::select(Assay, UniProt, cancer, sex, coef, se, FDR_BH) %>%
  dplyr::rename(P = FDR_BH,
         BETA = coef,
         SE = se) %>%
  dplyr::mutate(exposure_study = "EPIC",
         exposure_finemap = "InvRank")

df1 <- data_plot %>%
  dplyr::select(Assay, UniProt, cancer, sex, BETA, SE, P, exposure_study, exposure_finemap) 

data_plot_immuneonc <- bind_rows(df,df1) %>%
  dplyr::mutate(
    cancer = factor(cancer, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    exposure_finemap = factor(exposure_finemap, levels = c("InvRank", "coloc.abf", "p_ld-0.001"))
  ) %>%
  unique() %>%
  dplyr::group_by(UniProt, exposure_study, exposure_finemap, cancer, sex) %>%
  dplyr::mutate(row_index = case_when(
    exposure_study == "sun_2023_PMID37794186" ~ paste0("OlinkID: ", as.character(row_number())),
    exposure_study == "EPIC" ~ "EPIC",
    TRUE ~ NA_character_
  )) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Assay)  # Sort Assay alphabetically

# data_plot_explore ====
data_plot <- data_explore %>%
  dplyr::mutate(
    outcome = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    followup = factor(followup, levels = c(0, 2, 5)),
    coef_exp = as.numeric(coef_exp)
  ) %>%
  dplyr::arrange(exposure, outcome, sex, followup, model) %>%
  dplyr::group_by(exposure, outcome, sex, followup, model) %>%
  dplyr::mutate(FDR_BH = p.adjust(p = pval, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(followup == 0, model == "model_2") %>%
  dplyr::left_join(data_mr, by = c("UniProt", "outcome", "sex" = "outcome_sex")) %>%
  dplyr::filter(exposure_finemap %in% c("p_ld-0.001", "coloc.abf")) %>%
  dplyr::filter(outcome_study == "huyghe_2018_PMID30510241") %>%
  dplyr::filter(outcome_population == "ALL") %>%
  dplyr::filter(exposure_population == "ALL") %>%
  dplyr::mutate(
    cancer = factor(outcome, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male"))
  )

df <- data_plot %>%
  dplyr::select(Assay, UniProt, cancer, sex, coef, se, FDR_BH) %>%
  dplyr::rename(P = FDR_BH,
                BETA = coef,
                SE = se) %>%
  dplyr::mutate(exposure_study = "EPIC",
                exposure_finemap = "InvRank")

df1 <- data_plot %>%
  dplyr::select(Assay, UniProt, cancer, sex, BETA, SE, P, exposure_study, exposure_finemap) 

data_plot_explore <- bind_rows(df,df1) %>%
  dplyr::mutate(
    cancer = factor(cancer, levels = c("overall", "proximal", "distal")),
    sex = factor(sex, levels = c("combined", "female", "male")),
    exposure_finemap = factor(exposure_finemap, levels = c("InvRank", "coloc.abf", "p_ld-0.001"))
  ) %>%
  unique() %>%
  dplyr::group_by(UniProt, exposure_study, exposure_finemap, cancer, sex) %>%
  dplyr::mutate(row_index = case_when(
    exposure_study == "sun_2023_PMID37794186" ~ paste0("OlinkID: ", as.character(row_number())),
    exposure_study == "EPIC" ~ "EPIC",
    TRUE ~ NA_character_
  )) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Assay)  # Sort Assay alphabetically

# test data ====
temp <- data_plot_explore %>%
  dplyr::filter(Assay == "TFF3",
                cancer == "overall", 
                sex == "combined")

# plot: immuneonc ====
ID_assoc <- data_plot_immuneonc %>%
  dplyr::filter(exposure_study == "EPIC", P < 0.05) %>%
  dplyr::mutate(ID = paste0(UniProt, ";", cancer, ";", sex)) %>%
  dplyr::pull(ID)    

data_plot <- data_plot_immuneonc %>%
  dplyr::mutate(ID = paste0(UniProt, ";", cancer, ";", sex)) %>%
  dplyr::filter(ID %in% ID_assoc) 
  # dplyr::filter(
  #   (sex == "male" & cancer == "overall") |
  #     (sex == "female" & cancer == "overall") |
  #     (sex == "combined" & cancer %in% c("overall", "distal", "proximal"))
  # )

list_plots <- lapply(split(data_plot, list(data_plot$cancer, data_plot$sex), drop = TRUE), function(sub_data) {
  cancer_level <- unique(sub_data$cancer)
  sex_level <- unique(sub_data$sex)
  
  p <- functions::forestplot(df = sub_data, 
                             name = Assay, 
                             estimate = BETA, 
                             se = SE, 
                             pvalue = P, 
                             colour = exposure_study, 
                             shape = exposure_finemap, 
                             logodds = TRUE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    # ggplot2::coord_cartesian(xlim = c(min(sub_data$BETA), max(sub_data$BETA))) +
    ggplot2::labs(title = paste0("cancer = ", cancer_level, "; sex = ", sex_level),
                  colour = "study", shape = "method") +
    ggplot2::theme(legend.position = "right")
  
  # Save the plot
  outfile <- paste0("manuscript/figures/comparison_epic-mr/epic-immuneonc_", cancer_level, "_", sex_level, ".tiff")
  tiff(outfile, width = 600, height = 400, units = "px")
  print(p)
  dev.off()
  
  return(p)
})

# plot: explore ====
ID_assoc <- data_plot_explore %>%
  dplyr::filter(exposure_study == "EPIC", P < 0.05) %>%
  dplyr::mutate(ID = paste0(UniProt, ";", cancer, ";", sex)) %>%
  dplyr::pull(ID)    

data_plot <- data_plot_explore %>%
  dplyr::mutate(ID = paste0(UniProt, ";", cancer, ";", sex)) %>%
  dplyr::filter(ID %in% ID_assoc) 
# dplyr::filter(
#   (sex == "male" & cancer == "overall") |
#     (sex == "female" & cancer == "overall") |
#     (sex == "combined" & cancer %in% c("overall", "distal", "proximal"))
# )

list_plots <- lapply(split(data_plot, list(data_plot$cancer, data_plot$sex), drop = TRUE), function(sub_data) {
  cancer_level <- unique(sub_data$cancer)
  sex_level <- unique(sub_data$sex)
  
  p <- functions::forestplot(df = sub_data, 
                             name = Assay, 
                             estimate = BETA, 
                             se = SE, 
                             pvalue = P, 
                             colour = exposure_study, 
                             shape = exposure_finemap, 
                             logodds = TRUE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    ggplot2::coord_cartesian(xlim = c(min(exp(sub_data$BETA)), max(exp(sub_data$BETA)))) +
    ggplot2::labs(title = paste0("cancer = ", cancer_level, "; sex = ", sex_level),
                  colour = "study", shape = "method") +
    ggplot2::theme(legend.position = "right")
  
  # Save the plot
  outfile <- paste0("manuscript/figures/comparison_epic-mr/epic-explore_", cancer_level, "_", sex_level, ".tiff")
  tiff(outfile, width = 600, height = 1200, units = "px")
  print(p)
  dev.off()
  
  return(p)
})
