rm(list=ls())
set.seed(821)

# figure 2: forestplot of FDR associations for all followup levels

# environment ====
source("src/000_source.R")
rm(list=ls())

# data ====
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

# plot: immuneonc ====
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
  dplyr::filter(model == "model_2")

id_associations <- data_plot %>%
  dplyr::filter(
    followup == 0,
    model == "model_2",
    FDR_BH < 0.05
  ) %>%
  dplyr::select(exposure, outcome, sex) %>%
  dplyr::mutate(id = paste0(exposure, ";", outcome, ";", sex)) %>%
  dplyr::pull(id)

data_plot <- data_plot %>%
  dplyr::mutate(id = paste0(exposure, ";", outcome, ";", sex)) %>%
  dplyr::filter(id %in% id_associations)

## plot
combinations <- expand.grid(
  outcome = unique(data_plot$outcome),
  sex = unique(data_plot$sex),
  stringsAsFactors = FALSE
)

## plot
purrr::walk2(combinations$outcome, combinations$sex, function(cancer_level, sex_level) {
  
  # Filter data for this combination
  sub_data <- data_plot %>%
    dplyr::filter(outcome == cancer_level, sex == sex_level) %>%
    dplyr::arrange(Assay)
  
  # Generate the forest plot
  p <- functions::forestplot(df = sub_data, 
                             name = Assay, 
                             estimate = coef, 
                             se = se, 
                             pvalue = FDR_BH, 
                             colour = followup, 
                             logodds = TRUE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    ggplot2::labs(title = paste0(cancer_level, ";", sex_level),
                  colour = "study", shape = "method") +
    ggplot2::xlab("Odds ratio (95% confidence interval)") +
    ggplot2::theme(legend.position = "bottom")
  
  # Save the plot
  outfile <- paste0("manuscript/figures/followup/epic-immuneonc_", cancer_level, "_", sex_level, ".tiff")
  tiff(outfile, width = 1000, height = 1000, units = "px")
  print(p)
  dev.off()
})

# plot: explore ====
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
  dplyr::filter(model == "model_2")

id_associations <- data_plot %>%
  dplyr::filter(
    followup == 0,
    model == "model_2",
    FDR_BH < 0.05
  ) %>%
  dplyr::select(exposure, outcome, sex) %>%
  dplyr::mutate(id = paste0(exposure, ";", outcome, ";", sex)) %>%
  dplyr::pull(id)

data_plot <- data_plot %>%
  dplyr::mutate(id = paste0(exposure, ";", outcome, ";", sex)) %>%
  dplyr::filter(id %in% id_associations)

combinations <- expand.grid(
  outcome = unique(data_plot$outcome),
  sex = unique(data_plot$sex),
  stringsAsFactors = FALSE
)

## plot
purrr::walk2(combinations$outcome, combinations$sex, function(cancer_level, sex_level) {
  
  # Filter data for this combination
  sub_data <- data_plot %>%
    dplyr::filter(outcome == cancer_level, sex == sex_level) %>%
    dplyr::arrange(Assay)
  
  # Generate the forest plot
  p <- functions::forestplot(df = sub_data, 
                             name = Assay, 
                             estimate = coef, 
                             se = se, 
                             pvalue = FDR_BH, 
                             colour = followup, 
                             logodds = TRUE, 
                             psignif = 0.05, 
                             ci = 0.95) +
    ggplot2::labs(title = paste0(cancer_level, ";", sex_level),
         colour = "study", shape = "method") +
    ggplot2::xlab("Odds ratio (95% confidence interval)") +
    ggplot2::theme(legend.position = "bottom")
  
  # Save the plot
  outfile <- paste0("manuscript/figures/followup/epic-explore_", cancer_level, "_", sex_level, ".tiff")
  tiff(outfile, width = 1000, height = 1000, units = "px")
  print(p)
  dev.off()
})

