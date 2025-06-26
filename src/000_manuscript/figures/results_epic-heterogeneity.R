rm(list=ls())
set.seed(821)

# figure: volcano plot heterogeneity

# environment ====
source("src/000_source.R")
rm(list=ls())

# data: immuneonc ====
data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt")
data_info <- data.table::fread("data/processed/EPIC_olink-immuneonc.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% data$exposure)
data <- data %>%
  dplyr::left_join(data_info, by = c("exposure"="OlinkID")) %>%
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

data_sex <- data.table::fread("analysis/002_EPIC-analysis/002_heterogeneity/EPIC_olink-immuneonc_data-processed/heterogeneity-sex_data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt") %>%
  dplyr::filter(outcome == "overall", sex == "combined", model == "model_2", followup == 0) %>%
  dplyr::select(exposure, heterogeneity_p) %>%
  dplyr::rename(heterogeneity_sex_p = heterogeneity_p)

data_site <- data.table::fread("analysis/002_EPIC-analysis/002_heterogeneity/EPIC_olink-immuneonc_data-processed/heterogeneity-site_data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt") %>%
  dplyr::filter(outcome == "overall", sex == "combined", model == "model_2", followup == 0) %>%
  dplyr::select(exposure, heterogeneity_p) %>%
  dplyr::rename(heterogeneity_site_p = heterogeneity_p)

data_immuneonc <- data %>%
  dplyr::left_join(data_sex, by = "exposure") %>%
  dplyr::left_join(data_site, by = "exposure") %>%
  dplyr::mutate(
    heterogeneity = case_when(
      heterogeneity_sex_p < 0.05 & heterogeneity_site_p < 0.05 ~ "both",
      heterogeneity_sex_p < 0.05 ~ "sex",
      heterogeneity_site_p < 0.05 ~ "site",
      TRUE ~ NA_character_  # or "none" if you prefer
    )
  )

# data: explore ====
data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt")
data_info <- data.table::fread("data/processed/EPIC_olink-explore.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% data$exposure)
data <- data %>%
  dplyr::left_join(data_info, by = c("exposure"="OlinkID")) %>%
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

data_sex <- data.table::fread("analysis/002_EPIC-analysis/002_heterogeneity/EPIC_olink-explore_data-processed/heterogeneity-sex_data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt") %>%
  dplyr::filter(outcome == "overall", sex == "combined", model == "model_2", followup == 0) %>%
  dplyr::select(exposure, heterogeneity_p) %>%
  dplyr::rename(heterogeneity_sex_p = heterogeneity_p)

data_site <- data.table::fread("analysis/002_EPIC-analysis/002_heterogeneity/EPIC_olink-explore_data-processed/heterogeneity-site_data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt") %>%
  dplyr::filter(outcome == "overall", sex == "combined", model == "model_2", followup == 0) %>%
  dplyr::select(exposure, heterogeneity_p) %>%
  dplyr::rename(heterogeneity_site_p = heterogeneity_p)

data_explore <- data %>%
  dplyr::left_join(data_sex, by = "exposure") %>%
  dplyr::left_join(data_site, by = "exposure") %>%
  dplyr::mutate(
    heterogeneity = case_when(
      heterogeneity_sex_p < 0.05 & heterogeneity_site_p < 0.05 ~ "both",
      heterogeneity_sex_p < 0.05 ~ "sex",
      heterogeneity_site_p < 0.05 ~ "site",
      TRUE ~ NA_character_  # or "none" if you prefer
    )
  )

# plot: immuneonc - all ====
data_plot <- data_immuneonc
## plot
p1 <- ggplot2::ggplot() +
  
  # First set of points where heterogeneity is NA
  ggplot2::geom_point(
    data = data_plot %>%
      dplyr::filter(is.na(heterogeneity)), 
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH),
                 colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  
  # Second set of points where heterogeneity is NOT NA
  ggplot2::geom_point(
    data = data_plot %>%
      dplyr::filter(!is.na(heterogeneity)),  
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH),
                 fill = heterogeneity,colour = heterogeneity),
    shape = 21, size = 2) +
  
  # Lines
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") + 
  
  # Labels with correct dataset
  ggrepel::geom_text_repel(
    data = data_plot %>%
      dplyr::filter(!is.na(heterogeneity)),
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH), label = Assay),
    max.overlaps = 1000, min.segment.length = 0, force = 10,
    size = 5, box.padding = 0.5, point.padding = 0.5, colour = "black"
  ) +
  
  # Legend
  # FDR colour (first layer)
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("YES" = "orange", "NO" = "black", 
               "both" = "red3", "sex" = "lightgreen", "site" = "lightblue")) +
  # Fill colour for heterogeneity points
  ggplot2::scale_fill_manual(
    name = "Heterogeneity",
    values = c("both" = "red3", "sex" = "lightgreen", "site" = "lightblue")) +
  ggplot2::guides(alpha = "none", colour = "none", size = "none") +
  ggplot2::xlab("Odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  # Theme
  ggplot2::facet_grid(sex ~ outcome, scales = "fixed") +
  cowplot::theme_cowplot() 

tiff("manuscript/figures/heterogeneity/epic-immuneonc_all.tiff", 
     width = 1600, height = 800, units = "px")
p1
dev.off()

# plot: immuneonc - combined-overall ====
data_plot <- data_immuneonc %>%
  dplyr::filter(followup == 0, model == "model_2", sex == "combined", outcome == "overall")
## plot
p1 <- ggplot2::ggplot() +
  
  # First set of points where heterogeneity is NA
  ggplot2::geom_point(
    data = data_plot %>%
      dplyr::filter(is.na(heterogeneity)), 
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH),
                 colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  
  # Second set of points where heterogeneity is NOT NA
  ggplot2::geom_point(
    data = data_plot %>%
      dplyr::filter(!is.na(heterogeneity)),  
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH),
                 fill = heterogeneity,colour = heterogeneity),
    shape = 21, size = 2) +
  
  # Lines
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") + 
  
  # Labels with correct dataset
  ggrepel::geom_text_repel(
    data = data_plot %>%
      dplyr::filter(!is.na(heterogeneity)),
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH), label = Assay),
    max.overlaps = 1000, min.segment.length = 0, force = 10,
    size = 5, box.padding = 0.5, point.padding = 0.5, colour = "black"
  ) +
  
  # Legend
  # FDR colour (first layer)
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("YES" = "orange", "NO" = "black", 
               "both" = "red3", "sex" = "lightgreen", "site" = "lightblue")) +
  # Fill colour for heterogeneity points
  ggplot2::scale_fill_manual(
    name = "Heterogeneity",
    values = c("both" = "red3", "sex" = "lightgreen", "site" = "lightblue")) +
  ggplot2::guides(alpha = "none", colour = "none", size = "none") +
  ggplot2::xlab("Odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  # Theme
  ggplot2::facet_grid(sex ~ outcome, scales = "fixed") +
  cowplot::theme_cowplot() 

tiff("manuscript/figures/heterogeneity/epic-immuneonc_combined-overall.tiff", 
     width = 1600, height = 800, units = "px")
p1
dev.off()

# plot: explore - all ====
data_plot <- data_explore 
## plot
p1 <- ggplot2::ggplot() +
  
  # First set of points where heterogeneity is NA
  ggplot2::geom_point(
    data = data_plot %>%
      dplyr::filter(is.na(heterogeneity)), 
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH),
                 colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  
  # Second set of points where heterogeneity is NOT NA
  ggplot2::geom_point(
    data = data_plot %>%
      dplyr::filter(!is.na(heterogeneity)),  
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH),
                 fill = heterogeneity,colour = heterogeneity),
    shape = 21, size = 2) +
  
  # Lines
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") + 
  
  # Labels with correct dataset
  ggrepel::geom_text_repel(
    data = data_plot %>%
      dplyr::filter(!is.na(heterogeneity)),
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH), label = Assay),
    max.overlaps = 1000, min.segment.length = 0, force = 10,
    size = 5, box.padding = 0.5, point.padding = 0.5, colour = "black"
  ) +
  
  # Legend
  # FDR colour (first layer)
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("YES" = "orange", "NO" = "black", 
               "both" = "red3", "sex" = "lightgreen", "site" = "lightblue")) +
  # Fill colour for heterogeneity points
  ggplot2::scale_fill_manual(
    name = "Heterogeneity",
    values = c("both" = "red3", "sex" = "lightgreen", "site" = "lightblue")) +
  ggplot2::guides(alpha = "none", colour = "none", size = "none") +
  ggplot2::xlab("Odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  # Theme
  ggplot2::facet_grid(sex ~ outcome, scales = "fixed") +
  cowplot::theme_cowplot() 

## save
tiff("manuscript/figures/heterogeneity/epic-explore_all.tiff", 
     width = 1600, height = 800, units = "px")
p1
dev.off()

# plot: explore - combined-overall ====
data_plot <- data_explore %>%
  dplyr::filter(followup == 0, model == "model_2", sex == "combined", outcome == "overall")
## plot
p1 <- ggplot2::ggplot() +
  
  # First set of points where heterogeneity is NA
  ggplot2::geom_point(
    data = data_plot %>%
      dplyr::filter(is.na(heterogeneity)), 
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH),
                 colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  
  # Second set of points where heterogeneity is NOT NA
  ggplot2::geom_point(
    data = data_plot %>%
      dplyr::filter(!is.na(heterogeneity)),  
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH),
                 fill = heterogeneity,colour = heterogeneity),
    shape = 21, size = 2) +
  
  # Lines
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") + 
  
  # Labels with correct dataset
  ggrepel::geom_text_repel(
    data = data_plot %>%
      dplyr::filter(!is.na(heterogeneity)),
    ggplot2::aes(x = coef_exp, y = -log10(FDR_BH), label = Assay),
    max.overlaps = 1000, min.segment.length = 0, force = 10,
    size = 5, box.padding = 0.5, point.padding = 0.5, colour = "black"
  ) +
  
  # Legend
  # FDR colour (first layer)
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("YES" = "orange", "NO" = "black", 
               "both" = "red3", "sex" = "lightgreen", "site" = "lightblue")) +
  # Fill colour for heterogeneity points
  ggplot2::scale_fill_manual(
    name = "Heterogeneity",
    values = c("both" = "red3", "sex" = "lightgreen", "site" = "lightblue")) +
  ggplot2::guides(alpha = "none", colour = "none", size = "none") +
  ggplot2::xlab("Odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  # Theme
  ggplot2::facet_grid(sex ~ outcome, scales = "fixed") +
  cowplot::theme_cowplot() 

## save
tiff("manuscript/figures/heterogeneity/epic-explore_combined-overall.tiff", 
     width = 1600, height = 800, units = "px")
p1
dev.off()
