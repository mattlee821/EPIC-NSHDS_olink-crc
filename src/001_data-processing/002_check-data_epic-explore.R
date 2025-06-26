rm(list=ls())
set.seed(821)

# environment ====
source("src/000_source.R")
rm(list=ls())

# data ====
data <- data.table::fread("data/processed/EPIC_olink-explore.txt") %>%
  tibble::as_tibble() %>%
  dplyr::filter(SampleID != "Empty well")

df <- data %>%
  dplyr::select(SampleID, Match_Caseset, Sample_Type, PlateID, OlinkID, NPX) %>%
  unique() %>%
  tidyr::pivot_wider(names_from = OlinkID, values_from = NPX)

table_batch <- data %>%
  dplyr::select(SampleID, Sample_Type, PlateID, OlinkID, LOD) %>%
  unique() %>%
  dplyr::rename(feature = OlinkID)

table_all <- data %>%
  dplyr::select(SampleID, Sample_Type, PlateID, OlinkID, LOD) %>%
  unique() %>%
  dplyr::rename(feature = OlinkID)

# calculate missingness for SAMPLES and FEATURES ====
temp <- df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pct_missing_sample = mean(is.na(dplyr::c_across(dplyr::contains("OID")))) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::select(SampleID, pct_missing_sample)
table_batch <- table_batch %>%
  dplyr::left_join(temp, by = "SampleID")
table_all <- table_all %>%
  dplyr::left_join(temp, by = "SampleID")

temp <- df %>%
  dplyr::select(dplyr::contains("OID")) %>%
  sapply(function(col) mean(is.na(col)) * 100) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "feature") %>%
  tibble::as_tibble() %>%
  dplyr::rename(.data = ., pct_missing_feature = `.`)
table_batch <- table_batch %>%
  dplyr::left_join(temp, by = "feature")
table_all <- table_all %>%
  dplyr::left_join(temp, by = "feature")

# calculate W stat per FEATURE: for SAMPLES within batches ====
temp <- df %>%
  dplyr::filter(Sample_Type %in% c("SAMPLE", "CONTROL")) %>%
  dplyr::select(PlateID, Sample_Type, dplyr::contains("OID")) %>%
  tidyr::pivot_longer(
    cols = dplyr::contains("OID"),
    names_to = "feature",
    values_to = "value"
  ) %>%
  dplyr::group_by(PlateID, Sample_Type, feature) %>%
  dplyr::filter(dplyr::n() >= 3) %>%
  dplyr::summarise(
    W_stat = tryCatch({
      stats::shapiro.test(value)$statistic
    }, error = function(e) NA_real_),
    p_value = tryCatch({
      stats::shapiro.test(value)$p.value
    }, error = function(e) NA_real_),
    .groups = "drop"
  )
table_batch <- table_batch %>%
  dplyr::left_join(temp, by = c("PlateID", "Sample_Type", "feature"))

# calculate W stat per FEATURE: for SAMPLES across batches ====
temp <- df %>%
  dplyr::filter(Sample_Type %in% c("SAMPLE", "CONTROL")) %>%
  dplyr::select(Sample_Type, dplyr::contains("OID")) %>%
  tidyr::pivot_longer(
    cols = dplyr::contains("OID"),
    names_to = "feature",
    values_to = "value"
  ) %>%
  dplyr::group_by(Sample_Type, feature) %>%
  dplyr::filter(dplyr::n() >= 3) %>%
  dplyr::summarise(
    W_stat = tryCatch({
      stats::shapiro.test(value)$statistic
    }, error = function(e) NA_real_),
    p_value = tryCatch({
      stats::shapiro.test(value)$p.value
    }, error = function(e) NA_real_),
    .groups = "drop"
  )
table_all <- table_all %>%
  dplyr::left_join(temp, by = c("Sample_Type", "feature"))

# calculate CV for per FEATURE: within batches ====
temp <- df %>%
  dplyr::filter(Sample_Type %in% c("SAMPLE", "CONTROL")) %>%
  dplyr::select(PlateID, Sample_Type, dplyr::contains("OID")) %>%
  tidyr::pivot_longer(
    cols = dplyr::contains("OID"),
    names_to = "feature",
    values_to = "value"
  ) %>%
  dplyr::group_by(PlateID, Sample_Type, feature) %>%
  dplyr::summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = stats::sd(value, na.rm = TRUE),
    cv = ifelse(mean_val != 0, (sd_val / mean_val) * 100, NA_real_),
    n_samples = dplyr::n(),
    .groups = "drop"
  )
table_batch <- table_batch %>%
  dplyr::left_join(temp, by = c("PlateID", "Sample_Type", "feature"))

# calculate CV for per FEATURE: across batches ====
temp <- df %>%
  dplyr::filter(Sample_Type %in% c("SAMPLE", "CONTROL")) %>%
  dplyr::select(PlateID, Sample_Type, dplyr::contains("OID")) %>%
  tidyr::pivot_longer(
    cols = dplyr::contains("OID"),
    names_to = "feature",
    values_to = "value"
  ) %>%
  dplyr::group_by(Sample_Type, feature) %>%
  dplyr::summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = stats::sd(value, na.rm = TRUE),
    cv = ifelse(mean_val != 0, (sd_val / mean_val) * 100, NA_real_),
    n_samples = dplyr::n(),
    .groups = "drop"
  )
table_all <- table_all %>%
  dplyr::left_join(temp, by = c("Sample_Type", "feature"))

# ICC for per FEATURE: within batches - https://pmc.ncbi.nlm.nih.gov/articles/PMC6570933/ ====
temp <- df %>%
  dplyr::group_by(PlateID) %>%
  dplyr::group_map(.f = function(df_batch, batch_group) {
    batch <- batch_group$PlateID[1]
    
    # Extract sample and control data
    df_samples <- df_batch %>% dplyr::filter(Sample_Type == "SAMPLE")
    df_controls <- df_batch %>% dplyr::filter(Sample_Type == "CONTROL")
    
    # Skip if not enough data
    if (nrow(df_samples) < 2 || nrow(df_controls) < 2) return(NULL)
    
    # Pivot long, add index
    df_samples_long <- df_samples %>%
      dplyr::select(SampleID, Match_Caseset, dplyr::starts_with("OID")) %>%
      tidyr::pivot_longer(cols = dplyr::starts_with("OID"), names_to = "feature", values_to = "value") %>%
      dplyr::mutate(index = paste(SampleID, Match_Caseset, sep = "_")) %>%
      dplyr::select(feature, index, value) %>%
      tidyr::pivot_wider(names_from = index, values_from = value)
    
    df_controls_long <- df_controls %>%
      dplyr::select(SampleID, Match_Caseset, dplyr::starts_with("OID")) %>%
      tidyr::pivot_longer(cols = dplyr::starts_with("OID"), names_to = "feature", values_to = "value") %>%
      dplyr::mutate(index = paste(SampleID, Match_Caseset, sep = "_")) %>%
      dplyr::select(feature, index, value) %>%
      tidyr::pivot_wider(names_from = index, values_from = value)
    
    # Align rows on feature, then bind sample + control columns
    df_icc <- dplyr::inner_join(df_samples_long, df_controls_long, by = "feature") %>%
      tibble::column_to_rownames("feature")
    
    # Get final column names for ICC
    sample_cols <- colnames(df_samples_long)[-1]  # drop "feature"
    control_cols <- colnames(df_controls_long)[-1]
    
    # Run ICC
    icc_vals <- tryCatch(
      OmicsProcessing::calculate_ICC(
        df = df_icc,
        id_samples = sample_cols,
        id_qc = control_cols
      ),
      error = function(e) rep(NA_real_, nrow(df_icc))
    )
    
    tibble::tibble(
      feature = rownames(df_icc),
      PlateID = batch,
      ICC = icc_vals
    )
  }) %>%
  dplyr::bind_rows()
table_batch <- table_batch %>%
  dplyr::left_join(temp, by = c("PlateID", "feature"))

# ICC for per FEATURE: across batches - https://pmc.ncbi.nlm.nih.gov/articles/PMC6570933/ ====
temp <- {
  df_icc <- df %>%
    dplyr::select(SampleID, dplyr::contains("OID")) %>%
    t() %>%
    as.data.frame() %>%
    {\(x) {
      colnames(x) <- x[1, ]
      x[-1, , drop = FALSE]
    }}()
  
  id_samples <- df %>%
    dplyr::filter(Sample_Type == "SAMPLE") %>%
    dplyr::pull(SampleID) %>%
    unique()
  
  id_qc <- df %>%
    dplyr::filter(Sample_Type == "CONTROL") %>%
    dplyr::pull(SampleID) %>%
    unique()
  
  icc_vals <- tryCatch(
    OmicsProcessing::calculate_ICC(df = df_icc, id_samples = id_samples, id_qc = id_qc),
    error = function(e) rep(NA_real_, nrow(df_icc))
  )
  
  tibble::tibble(
    feature = rownames(df_icc),
    ICC = icc_vals
  )
}
table_all <- table_all %>%
  dplyr::left_join(temp, by = c("feature"))

# plot: distribution within batches ====
VAR_colour <- scales::hue_pal()(2)  # get first two default colors
grDevices::pdf("analysis/001_data-processing/EPIC_olink-explore_distribution-samples.pdf",
               width = 20, height = 20)

# Loop through each protein (OlinkID)
for (protein in unique(data$OlinkID)) {
  data_sub <- data %>% 
    dplyr::filter(OlinkID == protein)
  
  # skip if data_sub is empty or there are no SAMPLE entries
  if (nrow(data_sub) == 0 || !"SAMPLE" %in% data_sub$Sample_Type) {
    next
  }
  
  # combined plot
  p_combined <- ggplot2::ggplot(data_sub, 
                                ggplot2::aes(x = NPX, fill = Sample_Type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ PlateID, scales = "free") +
    cowplot::theme_cowplot() +
    ggplot2::labs(y = "Density", fill = element_blank()) +
    ggplot2::theme(legend.position = "top",
                   axis.title.x = element_blank())
  
  # samples plot
  p_samples <- ggplot2::ggplot(data_sub %>% 
                                 dplyr::filter(Sample_Type == "SAMPLE"),
                               ggplot2::aes(x = NPX, fill = Sample_Type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ PlateID, scales = "free") +
    cowplot::theme_cowplot() +
    ggplot2::labs(x = "NPX",y = "Density") +
    ggplot2::scale_fill_manual(values = VAR_colour[2]) +
    ggplot2::theme(legend.position = "none")
  
  # title plot
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      paste("Distribution across batches for: Assay", unique(data_sub$Assay), "; UniProt: ", unique(data_sub$UniProt), "; OlinkID:", protein),
      fontface = 'bold', x = 0, hjust = 0, size = 40) +
    ggplot2::theme(plot.margin = margin(0, 0, 0, 7))
  
  # plot
  p <- cowplot::plot_grid(title, p_combined, p_samples, 
                          ncol = 1, rel_heights = c(0.1,1,1))
  print(p)
}

# Close the PDF device
grDevices::dev.off()

# plot: distribution across batches  ====
# Order IDs by PlateID and SampleID
id_order <- data %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::arrange(PlateID, SampleID) %>%
  dplyr::pull(SampleID) %>%
  unique()

data_plot <- data %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::mutate(SampleID = factor(SampleID, levels = id_order))

grDevices::pdf("analysis/001_data-processing/EPIC_olink-explore_scatter-batch.pdf",
               width = 16, height = 8)

# Loop through each protein (OlinkID)
for (protein in unique(data_plot$OlinkID)) {
  
  data_sub <- data_plot %>%
    dplyr::filter(OlinkID == protein)
  
  # skip if data_sub is empty or there are no SAMPLE entries
  if (nrow(data_sub) == 0 || !"SAMPLE" %in% data_sub$Sample_Type) {
    next
  }
  
  # 1. Create a batch-level summary dataframe using external info
  batch_summary <- table_batch %>%
    dplyr::filter(feature == protein & Sample_Type == "SAMPLE") %>%
    dplyr::group_by(PlateID) %>%
    dplyr::summarise(
      missing_feature = unique(pct_missing_feature),
      LOD = unique(LOD),
      CV = unique(cv),
      ICC = unique(ICC)
    )
  
  batches <- unique(data_sub$PlateID)
  batch_colors <- scales::hue_pal()(length(batches))
  names(batch_colors) <- batches
  
  make_colored_line <- function(prefix, NPXs) {
    purrr::map2_chr(NPXs, batch_colors[batches], ~ {
      val <- ifelse(is.na(.x), "NA", format(round(.x, 2), nsmall = 2))
      paste0("<span style='color:", .y, "'>", val, "</span>")
    }) %>%
      paste(collapse = "; ") %>%
      paste0("<b>", prefix, ":</b> ", .)
  }
  
  missing_line <- make_colored_line("Missing", batch_summary$missing_feature)
  LOD_line <- make_colored_line("LOD", batch_summary$LOD)
  CV_line <- make_colored_line("CV", batch_summary$CV)
  ICC_line <- make_colored_line("ICC", batch_summary$ICC)
  
  full_label <- paste(missing_line, 
                      LOD_line,
                      CV_line,
                      ICC_line,
                      sep = "<br>")
  
  # scatter plot
  p <- ggplot2::ggplot(data_sub, 
                       ggplot2::aes(
                         x = SampleID, y = NPX, 
                         color = PlateID, shape = as.factor(Cncr_Caco_Clrt))) +
    ggplot2::geom_point(alpha = 1) +
    ggplot2::labs(
      title = paste("Sample NPX values across batches for: Assay", unique(data_sub$Assay), "; UniProt: ", unique(data_sub$UniProt), "; OlinkID:", protein),
      color = "Batch", shape = "Cancer"
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "vertical",  # stack legends vertically
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = ceiling(length(batches)/5))) +
    ggplot2::geom_hline(data = batch_summary,
                        ggplot2::aes(yintercept = LOD, color = PlateID),
                        linewidth = 1.5, linetype = "solid")
  
  # label plot
  label_plot <- ggplot2::ggplot() +
    ggtext::geom_richtext(
      data = data.frame(x = 0.03, y = 0.8, label = full_label),
      ggplot2::aes(x = x, y = y, label = label),
      hjust = 0, vjust = 1,
      size = 5,
      label.color = NA,
      fill = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1), expand = FALSE, clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"))
  
  # plot
  p <- cowplot::plot_grid(
    p,label_plot,
    ncol = 1, rel_heights = c(1,0.2)
  )
  
  print(p)
}

# Close the PDF device
grDevices::dev.off()

# write ====
write.table(table_all, "analysis/001_data-processing/EPIC_olink-explore_metadata-across-batches.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(table_batch, "analysis/001_data-processing/EPIC_olink-explore_metadata-within-batches.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
