rm(list=ls())
set.seed(821)

# section 1: results figures

# environment ====
source("src/000_source.R")
rm(list=ls())

# data ====
data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_intersect_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-FALSE_transformation-InvRank_outlier-FALSE_platecorrection-FALSE_centre-scale-TRUE.txt")
data_info <- data.table::fread("data/processed/EPIC_olink-explore.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% data$exposure)
data_explore <- data %>%
  dplyr::left_join(data_info, by = c("exposure"="OlinkID")) %>%
  tibble::as_tibble() %>%
  dplyr::group_by(UniProt) %>%
  dplyr::mutate(exposure_index = as.factor(match(exposure, unique(exposure)))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(analysis = "explore",
                coef_exp = as.numeric(coef_exp),
                ci_lower_exp = as.numeric(ci_lower_exp),
                ci_upper_exp = as.numeric(ci_upper_exp)) %>%
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
  dplyr::filter(followup == 0, model == "model_2")

data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_intersect_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-FALSE_transformation-InvRank_outlier-FALSE_platecorrection-FALSE_centre-scale-TRUE.txt")
data_info <- data.table::fread("data/processed/EPIC_olink-immuneonc.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% data$exposure)
data_immuneonc <- data %>%
  dplyr::left_join(data_info, by = c("exposure"="OlinkID")) %>%
  tibble::as_tibble() %>%
  dplyr::group_by(UniProt) %>%
  dplyr::mutate(exposure_index = as.factor(match(exposure, unique(exposure)))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(analysis = "immuneonc",
                coef_exp = as.numeric(coef_exp),
                ci_lower_exp = as.numeric(ci_lower_exp),
                ci_upper_exp = as.numeric(ci_upper_exp)) %>%
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
  dplyr::filter(followup == 0, model == "model_2")

ID_immuneonc <- unique(data_immuneonc$UniProt)
ID_explore <- unique(data_immuneonc$UniProt)

# plot all ====
data_plot <- dplyr::bind_rows(data_explore, data_immuneonc)  %>%
  dplyr::group_by(UniProt) %>% # match the Assay value for explore and immuneonc using immuneonc value
  dplyr::mutate(
    Assay = if_else(
      analysis == "explore",
      Assay[analysis == "immuneonc"][1],  # get immuneonc Assay for this UniProt
      Assay
    )
  ) %>%
  dplyr::ungroup()

list_plots <- lapply(split(data_plot, list(data_plot$outcome, data_plot$sex), drop = TRUE), function(sub_data) {
  cancer_level <- unique(sub_data$outcome)
  sex_level <- unique(sub_data$sex)
  
  p <- functions::forestplot(df = sub_data, 
                             name = Assay, 
                             estimate = coef, 
                             se = se, 
                             pvalue = FDR_BH, 
                             colour = analysis, 
                             shape = exposure_index, 
                             logodds = TRUE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    ggplot2::coord_cartesian(xlim = c(min(exp(sub_data$coef)), max(exp(sub_data$coef)))) +
    ggplot2::labs(title = paste0("cancer = ", cancer_level, "; sex = ", sex_level),
                  colour = "study", shape = "index") +
    ggplot2::theme(legend.position = "right")
  
  # Save the plot
  outfile <- paste0("manuscript/figures/comparison_epic-intersect/", cancer_level, "_", sex_level, ".tiff")
  tiff(outfile, width = 600, height = 1200, units = "px")
  print(p)
  dev.off()
  
  return(p)
})

# plot sig all ====
ID_explore <- data_explore %>%
  dplyr::filter(FDR_BH < 0.05) %>%
  dplyr::mutate(ID = paste0(UniProt, ";", outcome, ";", sex)) %>%
  dplyr::pull(ID)    
ID_immuneonc <- data_immuneonc %>%
  dplyr::filter(FDR_BH < 0.05) %>%
  dplyr::mutate(ID = paste0(UniProt, ";", outcome, ";", sex)) %>%
  dplyr::pull(ID)   
ID_all <- c(unique(ID_explore), unique(ID_immuneonc))

data_plot <- dplyr::bind_rows(data_explore, data_immuneonc)  %>%
  dplyr::group_by(UniProt) %>% # match the Assay value for explore and immuneonc using immuneonc value
  dplyr::mutate(
    Assay = if_else(
      analysis == "explore",
      Assay[analysis == "immuneonc"][1],  # get immuneonc Assay for this UniProt
      Assay
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ID = paste0(UniProt, ";", outcome, ";", sex)) %>%
  dplyr::filter(ID %in% ID_all) 

list_plots <- lapply(split(data_plot, list(data_plot$outcome, data_plot$sex), drop = TRUE), function(sub_data) {
  cancer_level <- unique(sub_data$outcome)
  sex_level <- unique(sub_data$sex)
  
  p <- functions::forestplot(df = sub_data, 
                             name = Assay, 
                             estimate = coef, 
                             se = se, 
                             pvalue = FDR_BH, 
                             colour = analysis, 
                             shape = exposure_index, 
                             logodds = TRUE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    # ggplot2::coord_cartesian(xlim = c(min(exp(sub_data$coef)), max(exp(sub_data$coef)))) +
    ggplot2::labs(title = paste0("cancer = ", cancer_level, "; sex = ", sex_level),
                  colour = "study", shape = "index") +
    ggplot2::theme(legend.position = "right")
  
  # Save the plot
  outfile <- paste0("manuscript/figures/comparison_epic-intersect/sig-all_", cancer_level, "_", sex_level, ".tiff")
  tiff(outfile, width = 600, height = 1200, units = "px")
  print(p)
  dev.off()
  
  return(p)
})

# plot sig shared ====
ID_explore <- data_explore %>%
  dplyr::filter(FDR_BH < 0.05) %>%
  dplyr::mutate(ID = paste0(UniProt, ";", outcome, ";", sex)) %>%
  dplyr::pull(ID)    
ID_immuneonc <- data_immuneonc %>%
  dplyr::filter(FDR_BH < 0.05) %>%
  dplyr::mutate(ID = paste0(UniProt, ";", outcome, ";", sex)) %>%
  dplyr::pull(ID)   
ID_shared <- intersect(ID_explore, ID_immuneonc)

data_plot <- dplyr::bind_rows(data_explore, data_immuneonc)  %>%
  dplyr::group_by(UniProt) %>% # match the Assay value for explore and immuneonc using immuneonc value
  dplyr::mutate(
    Assay = if_else(
      analysis == "explore",
      Assay[analysis == "immuneonc"][1],  # get immuneonc Assay for this UniProt
      Assay
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ID = paste0(UniProt, ";", outcome, ";", sex)) %>%
  dplyr::filter(ID %in% ID_shared) 

list_plots <- lapply(split(data_plot, list(data_plot$outcome, data_plot$sex), drop = TRUE), function(sub_data) {
  cancer_level <- unique(sub_data$outcome)
  sex_level <- unique(sub_data$sex)
  
  p <- functions::forestplot(df = sub_data, 
                             name = Assay, 
                             estimate = coef, 
                             se = se, 
                             pvalue = FDR_BH, 
                             colour = analysis, 
                             shape = exposure_index, 
                             logodds = TRUE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    # ggplot2::coord_cartesian(xlim = c(min(exp(sub_data$coef)), max(exp(sub_data$coef)))) +
    ggplot2::labs(title = paste0("cancer = ", cancer_level, "; sex = ", sex_level),
                  colour = "study", shape = "index") +
    ggplot2::theme(legend.position = "right")
  
  # Save the plot
  outfile <- paste0("manuscript/figures/comparison_epic-intersect/sig-shared_", cancer_level, "_", sex_level, ".tiff")
  tiff(outfile, width = 600, height = 1200, units = "px")
  print(p)
  dev.off()
  
  return(p)
})
