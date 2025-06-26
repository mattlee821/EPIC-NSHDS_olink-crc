rm(list=ls())
set.seed(821)

# figure 1: volcano plot all

# environment ====
source("src/000_source.R")
rm(list=ls())

# data ====
data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-explore_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt")
data_info_explore <- data.table::fread("data/processed/EPIC_olink-explore.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% data$exposure)
data_explore <- data %>%
  dplyr::left_join(data_info_explore, by = c("exposure"="OlinkID"))

data <- data.table::fread("analysis/002_EPIC-analysis/001_analysis/EPIC_olink-immuneonc_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-LOD_transformation-InvRank_outlier-TRUE_platecorrection-FALSE_centre-scale-TRUE.txt")
data_info_immuneonc <- data.table::fread("data/processed/EPIC_olink-immuneonc.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% data$exposure)
data_immuneonc <- data %>%
  dplyr::left_join(data_info_immuneonc, by = c("exposure"="OlinkID"))

# mapping ====
data_info <- data_info_explore %>%
  dplyr::select(UniProt, Assay, OlinkID) %>%
  unique() %>%
  dplyr::left_join(data_info_immuneonc %>%
                     dplyr::select(UniProt, Assay, OlinkID) %>%
              unique(),
            by = "UniProt") %>%
  dplyr::rename(Assay_explore = Assay.x,
         Assay_immuneonc = Assay.y,
         OlinkID_explore = OlinkID.x,
         OlinkID_immuneonc = OlinkID.y) %>%
  # Create unified OlinkID column for indexing
  dplyr::mutate(OlinkID = coalesce(OlinkID_explore, OlinkID_immuneonc),
         Assay_base = coalesce(Assay_explore, Assay_immuneonc)) %>%
  dplyr::group_by(Assay_base) %>%
  dplyr::mutate(index = row_number(),
         Assay_map = if (n() > 1) paste0(Assay_base, "_", index) else Assay_base) %>%
  dplyr::ungroup() %>%
  dplyr::select(-OlinkID, -Assay_base, -index)

# format ====
data_immuneonc <- data_immuneonc %>%
  dplyr::left_join(data_info %>%
                     dplyr::select(UniProt, OlinkID_immuneonc, Assay_explore) %>%
                     unique(), 
                   by = c("exposure" = "OlinkID_immuneonc", "UniProt" = "UniProt"))
data_explore <- data_explore %>%
  dplyr::left_join(data_info %>%
                     dplyr::select(UniProt, OlinkID_explore, Assay_map) %>%
                     unique(), 
                   by = c("exposure" = "OlinkID_explore", "UniProt" = "UniProt"))

# plot: immuneonc - all ====
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
  dplyr::filter(followup == 0, model == "model_2")

p1 <- ggplot2::ggplot(data_plot, 
                      ggplot2::aes(x = coef_exp, y = -log10(FDR_BH))) +
  ggplot2::geom_point(
    ggplot2::aes(colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  ggplot2::scale_color_manual(values = c("YES" = "orange", "NO" = "black")) +
  cowplot::theme_cowplot() +
  ggplot2::facet_grid(sex~outcome, scales = "fixed") +
  
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") +
  
  # legend
  ggplot2::guides(alpha = "none", size = "none", colour = "none") +
  # labels
  ggplot2::xlab("odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +

  ggrepel::geom_text_repel(
    data = data_plot %>%
      filter(FDR_BH < 0.05),
    ggplot2::aes(x = coef_exp,
                 y = -log10(FDR_BH),
                 label = Assay_explore),
    max.overlaps = 1000,
    min.segment.length = 0,
    colour = "black",
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5
)

## save
tiff("manuscript/figures/results/epic-immuneonc_all.tiff", 
     width = 1400, height = 600, units = "px")
p1
dev.off()

# plot: immuneonc - combined ====
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
  dplyr::filter(followup == 0, model == "model_2", sex == "combined")

p1 <- ggplot2::ggplot(data_plot, 
                      ggplot2::aes(x = coef_exp, y = -log10(FDR_BH))) +
  ggplot2::geom_point(
    ggplot2::aes(colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  ggplot2::scale_color_manual(values = c("YES" = "orange", "NO" = "black")) +
  cowplot::theme_cowplot() +
  ggplot2::facet_grid(~outcome, scales = "fixed") +
  
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") +
  
  # legend
  ggplot2::guides(alpha = "none", size = "none", colour = "none") +
  # labels
  ggplot2::xlab("odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  ggrepel::geom_text_repel(
    data = data_plot %>%
      filter(FDR_BH < 0.05),
    ggplot2::aes(x = coef_exp,
                 y = -log10(FDR_BH),
                 label = Assay_explore),
    max.overlaps = 1000,
    min.segment.length = 0,
    colour = "black",
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5
  )

## save
tiff("manuscript/figures/results/epic-immuneonc_combined.tiff", 
     width = 1400, height = 400, units = "px")
p1
dev.off()

# plot: immuneonc - combined-overall ====
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
  dplyr::filter(followup == 0, model == "model_2", sex == "combined", outcome == "overall")

p1 <- ggplot2::ggplot(data_plot, 
                      ggplot2::aes(x = coef_exp, y = -log10(FDR_BH))) +
  ggplot2::geom_point(
    ggplot2::aes(colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  ggplot2::scale_color_manual(values = c("YES" = "orange", "NO" = "black")) +
  cowplot::theme_cowplot() +

  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") +
  
  # legend
  ggplot2::guides(alpha = "none", size = "none", colour = "none") +
  # labels
  ggplot2::xlab("odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  ggrepel::geom_text_repel(
    data = data_plot %>%
      filter(FDR_BH < 0.05),
    ggplot2::aes(x = coef_exp,
                 y = -log10(FDR_BH),
                 label = Assay_explore),
    max.overlaps = 1000,
    min.segment.length = 0,
    colour = "black",
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5
  )

## save
tiff("manuscript/figures/results/epic-immuneonc_combined-overall.tiff", 
     width = 800, height = 400, units = "px")
p1
dev.off()

# plot: explore - all ====
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
  dplyr::filter(followup == 0, model == "model_2")

p1 <- ggplot2::ggplot(data_plot, 
                      ggplot2::aes(x = coef_exp, y = -log10(FDR_BH))) +
  ggplot2::geom_point(
    ggplot2::aes(colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  ggplot2::scale_color_manual(values = c("YES" = "orange", "NO" = "black")) +
  cowplot::theme_cowplot() +
  ggplot2::facet_grid(sex~outcome, scales = "fixed") +
  
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") +
  
  # legend
  ggplot2::guides(alpha = "none", size = "none", colour = "none") +
  # labels
  ggplot2::xlab("odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  ggrepel::geom_text_repel(
    data = data_plot %>%
      filter(FDR_BH < 0.05),
    ggplot2::aes(x = coef_exp,
                 y = -log10(FDR_BH),
                 label = Assay_map),
    max.overlaps = 1000,
    min.segment.length = 0,
    colour = "black",
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5
  )

## save
tiff("manuscript/figures/results/epic-explore_all.tiff", 
     width = 1400, height = 600, units = "px")
p1
dev.off()

# plot: explore - combined ====
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
  dplyr::filter(followup == 0, model == "model_2", sex == "combined")

p1 <- ggplot2::ggplot(data_plot, 
                      ggplot2::aes(x = coef_exp, y = -log10(FDR_BH))) +
  ggplot2::geom_point(
    ggplot2::aes(colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  ggplot2::scale_color_manual(values = c("YES" = "orange", "NO" = "black")) +
  cowplot::theme_cowplot() +
  ggplot2::facet_grid(~outcome, scales = "fixed") +
  
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") +
  
  # legend
  ggplot2::guides(alpha = "none", size = "none", colour = "none") +
  # labels
  ggplot2::xlab("odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  ggrepel::geom_text_repel(
    data = data_plot %>%
      filter(FDR_BH < 0.05),
    ggplot2::aes(x = coef_exp,
                 y = -log10(FDR_BH),
                 label = Assay_map),
    max.overlaps = 1000,
    min.segment.length = 0,
    colour = "black",
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5
  )

## save
tiff("manuscript/figures/results/epic-explore_combined.tiff", 
     width = 1600, height = 400, units = "px")
p1
dev.off()

# plot: explore - combined-overall ====
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
  dplyr::filter(followup == 0, model == "model_2", sex == "combined", outcome == "overall")

p1 <- ggplot2::ggplot(data_plot, 
                      ggplot2::aes(x = coef_exp, y = -log10(FDR_BH))) +
  ggplot2::geom_point(
    ggplot2::aes(colour = ifelse(FDR_BH < 0.05, "YES", "NO"))) +
  ggplot2::scale_color_manual(values = c("YES" = "orange", "NO" = "black")) +
  cowplot::theme_cowplot() +
  
  ggplot2::geom_vline(xintercept = 1, color = "grey", linetype = "solid") +
  ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") +
  
  # legend
  ggplot2::guides(alpha = "none", size = "none", colour = "none") +
  # labels
  ggplot2::xlab("odds ratio") +
  ggplot2::ylab(label = "-log10(FDR pval)") +
  
  ggrepel::geom_text_repel(
    data = data_plot %>%
      filter(FDR_BH < 0.05),
    ggplot2::aes(x = coef_exp,
                 y = -log10(FDR_BH),
                 label = Assay_map),
    max.overlaps = 1000,
    min.segment.length = 0,
    colour = "black",
    size = 5,
    box.padding = 0.5,
    point.padding = 0.5
  )

## save
tiff("manuscript/figures/results/epic-explore_combined-overall.tiff", 
     width = 800, height = 400, units = "px")
p1
dev.off()
