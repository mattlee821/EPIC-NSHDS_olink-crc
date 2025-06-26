rm(list=ls())
set.seed(821)

#environment ====
source("src/000_source.R")
rm(list=ls())

# data ====
load("analysis/001_data-processing/EPIC_olink-immuneonc_intersect_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-FALSE_transformation-InvRank_outlier-FALSE_platecorrection-FALSE_centre-scale-TRUE.rda")
data_immuneonc <- df_processed %>%
  dplyr::arrange(SampleID) %>%
  tibble::column_to_rownames(var = "SampleID")

load("analysis/001_data-processing/EPIC_olink-explore_intersect_data-processed/data-features_exclusion-feature-0.8_exclusion-sample-0.8_imputation-FALSE_transformation-InvRank_outlier-FALSE_platecorrection-FALSE_centre-scale-TRUE.rda")
data_explore <- df_processed %>%
  dplyr::arrange(SampleID) %>%
  tibble::column_to_rownames(var = "SampleID")

stopifnot(identical(rownames(data_explore), rownames(data_immuneonc)))

# mapping ====
map_explore <- data.table::fread("data/processed/EPIC_olink-explore.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% colnames(data_explore))

map_immuneonc <- data.table::fread("data/processed/EPIC_olink-immuneonc.txt") %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
  unique() %>%
  dplyr::filter(OlinkID %in% colnames(data_immuneonc)) %>%
  distinct(UniProt, OlinkID) %>%
  dplyr::filter(UniProt %in% map_explore$UniProt)

# format: rename duplicated Assays in explore with indexed suffixes ====
map_explore <- map_explore %>%
  dplyr::group_by(UniProt) %>%
  dplyr::mutate(UniProt_indexed = if (n() == 1) UniProt else paste0(UniProt, "_", row_number())) %>%
  dplyr::ungroup()
# Rename columns in data_explore
colnames(data_explore)[match(map_explore$OlinkID, colnames(data_explore))] <- map_explore$UniProt_indexed
# Rename columns in data_immuneonc to Assay (simple names)
colnames(data_immuneonc)[match(map_immuneonc$OlinkID, colnames(data_immuneonc))] <- map_immuneonc$UniProt

# cor ====
uniptot_base_names <- gsub("_[0-9]+$", "", colnames(data_explore))
df_cor <- purrr::map2_dfr(
  .x = colnames(data_explore),
  .y = uniptot_base_names,
  .f = function(explore_col, immuneonc_col) {
    x <- data_explore[[explore_col]]
    y <- data_immuneonc[[immuneonc_col]]
    
    tibble(
      explore_col = explore_col,
      immuneonc_col = immuneonc_col,
      spearman = cor(x, y, method = "spearman", use = "complete.obs"),
      pearson = cor(x, y, method = "pearson", use = "complete.obs")
    ) %>%
      dplyr::left_join(map_explore %>%
                  dplyr::select(UniProt, Assay), 
                  by = c("immuneonc_col" = "UniProt")) %>%
      unique() %>%
      dplyr::select(Assay, immuneonc_col, spearman, pearson) %>%
      dplyr::rename(UniProt = immuneonc_col)
  }
)

# format df ====
df_cor <- dplyr::bind_rows(df_cor) %>%
  tidyr::pivot_longer(cols = c(pearson, spearman), names_to = "Method", values_to = "Correlation")
df_cor %>%
  group_by(Method) %>%
  summarise(
    mean = mean(Correlation, na.rm = TRUE),
    median = median(Correlation, na.rm = TRUE),
    sd = sd(Correlation, na.rm = TRUE),
    min = min(Correlation, na.rm = TRUE),
    max = max(Correlation, na.rm = TRUE),
    .groups = "drop"
  )

# plot: scatter for each protein ====
# Convert rownames to column "SampleID"
data_explore <- data_explore %>% rownames_to_column(var = "SampleID")
data_immuneonc <- data_immuneonc %>% rownames_to_column(var = "SampleID")
# Get base UniProt names (without index suffix) to match with immuneonc
uniptot_base_names <- gsub("_[0-9]+$", "", map_explore$UniProt_indexed)
# Create a dataframe to guide plotting
df_plot_map <- tibble(
  explore_col = map_explore$UniProt_indexed,
  immuneonc_col = uniptot_base_names,
  Assay = map_explore$Assay
) 
# Make plots
plot_list <- purrr::pmap(
  .l = list(
    explore_col = df_plot_map$explore_col,
    immuneonc_col = df_plot_map$immuneonc_col,
    Assay = df_plot_map$Assay
  ),
  .f = function(explore_col, immuneonc_col, Assay) {
    
    df_tmp <- data_explore %>%
      dplyr::select(SampleID, Explore = all_of(explore_col)) %>%
      dplyr::left_join(
        data_immuneonc %>%
          dplyr::select(SampleID, ImmuneOnc = all_of(immuneonc_col)),
        by = "SampleID"
      )
    
    # Correlations
    pearson_cor <- cor(df_tmp$Explore, df_tmp$ImmuneOnc, method = "pearson", use = "complete.obs")
    spearman_cor <- cor(df_tmp$Explore, df_tmp$ImmuneOnc, method = "spearman", use = "complete.obs")
    
    # plot
    ggplot2::ggplot(df_tmp, 
                    ggplot2::aes(x = ImmuneOnc, y = Explore)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 0, color = "grey70") +
      ggplot2::geom_hline(yintercept = 0, color = "grey70") +
      ggplot2::labs(
        title = paste0(Assay),
        x = paste0("immuneonc: ", immuneonc_col),
        y = paste0("explore: ", explore_col),
        subtitle = paste0(
          "Pearson: ", round(pearson_cor, 2),
          " | Spearman: ", round(spearman_cor, 2)
        )
      ) +
      cowplot::theme_cowplot()
  }
)

# Save all plots in one PDF
pdf("manuscript/figures/comparison_epic-intersect/scatterplots_explore_vs_immuneonc.pdf",
    width = 10, height = 10)

for (p in plot_list) {
  print(p)
}

dev.off()

# plot: density/box ====
# Compute summary statistics
summary_stats <- df_cor %>%
  group_by(Method) %>%
  summarise(
    mean = mean(Correlation, na.rm = TRUE),
    sd = sd(Correlation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(mean, sd), names_to = "stat", values_to = "value") %>%
  tidyr::pivot_wider(names_from = stat, values_from = value)

# Create plot
p1 <- ggplot2::ggplot(df_cor, ggplot2::aes(x = Correlation, fill = Method)) +
  ggplot2::geom_density(alpha = 0.5) +
  ggplot2::geom_vline(data = summary_stats, aes(xintercept = mean, color = Method), linetype = "solid") +
  ggplot2::geom_vline(data = summary_stats, aes(xintercept = mean + sd, color = Method), linetype = "dashed") +
  ggplot2::geom_vline(data = summary_stats, aes(xintercept = mean - sd, color = Method), linetype = "dashed") +
  cowplot::theme_cowplot() +
  ggplot2::labs(
    title = "Distribution of Correlations",
    x = "Correlation Coefficient", y = "Density"
  ) 

# plot: boxplot
p2 <- ggplot(df_cor, ggplot2::aes(x = Method, y = Correlation, fill = Method)) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  cowplot::theme_cowplot() +
  ggplot2::labs(
    title = "Boxplot of correlations",
    x = "Correlation Method",
    y = "Correlation Coefficient"
  ) +
  ggplot2::scale_y_continuous(limits = c(min(df_cor$Correlation)*1.1, 1)) +
  ggplot2::theme(legend.position = "none")

# write 
grDevices::tiff(filename = "manuscript/figures/comparison_epic-intersect/density-box_explore_vs_immuneonc.tiff", 
                width = 800, height = 600, units = "px")
cowplot::plot_grid(p1,p2, nrow = 1)
grDevices::dev.off()

