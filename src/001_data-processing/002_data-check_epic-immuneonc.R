rm(list=ls())
set.seed(821)

# environment ====

# data ====
data_meta_controls <- data.table::fread("data/processed/EPIC_olink-immuneonc_sample-metadata-controls.txt")
data_meta_features <- data.table::fread("data/processed/EPIC_olink-immuneonc_feature-metadata.txt")
data_meta_samples <- data.table::fread("data/processed/EPIC_olink-immuneonc_sample-metadata.txt")
data_features_controls <- data.table::fread("data/processed/EPIC_olink-immuneonc_feature-data-controls.txt")
data_features_samples <- data.table::fread("data/processed/EPIC_olink-immuneonc_feature-data.txt")

# calculate sample missing ====
data_meta_samples <- data_meta_samples %>%
  dplyr::left_join(
    data_features_samples %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        missing_pct = sum(is.na(dplyr::c_across(-Idepic_Bio))) / (ncol(data_features_samples) - 1)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(Idepic_Bio, missing_pct), 
    by = "Idepic_Bio"
  )
write.table(data_meta_samples, "analysis/001_data-processing/EPIC_olink-immuneonc_sample-metadata.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# W statistics of normality ====
data_meta_features$W_stat_samples <- apply(data_features_samples[, -1], 2, function(protein) {
  if (all(is.na(protein))) {
    return(NA)
  }
  if (length(unique(protein[!is.na(protein)])) > 3) {
    stats::shapiro.test(protein)$statistic
  } else {
    NA
  }
})

data_meta_features$W_stat_controls <- apply(data_features_controls[, -1], 2, function(protein) {
  if (all(is.na(protein))) {
    return(NA)
  }
  if (length(unique(protein[!is.na(protein)])) > 3) {
    stats::shapiro.test(protein)$statistic
  } else {
    NA
  }
})

# calculate CV for each protein ====
data_meta_features <- data_meta_features %>%
  dplyr::left_join(
    data_features_samples %>%
      dplyr::select(-Idepic_Bio) %>%  # Exclude the Idepic_Bio column
      dplyr::summarise(dplyr::across(everything(), ~ {
        sd_value <- sd(.x, na.rm = TRUE)  # Calculate standard deviation, excluding NAs
        mean_value <- mean(.x, na.rm = TRUE)  # Calculate mean, excluding NAs
        cv_value <- (sd_value / mean_value) * 100  # Calculate CV
        return(cv_value)  # Return the CV value
      })) %>%
      tidyr::pivot_longer(everything(), names_to = "OlinkID", values_to = "CV"), 
    by = "OlinkID"  # Join by UNIPROT
  )

# calculate ICC for each protein - https://pmc.ncbi.nlm.nih.gov/articles/PMC6570933/ ====
id_samples <- data_features_samples$Idepic_Bio
id_qc <- data_features_controls$Idepic_Bio
data_icc <- data_features_samples %>%
  tidyr::pivot_longer(cols = -Idepic_Bio, names_to = "OlinkID", values_to = "value") %>%  # Convert samples to long format
  tidyr::pivot_wider(names_from = Idepic_Bio, values_from = value) %>%  # Transpose samples
  dplyr::left_join(
    data_features_controls %>%
      tidyr::pivot_longer(cols = -Idepic_Bio, names_to = "OlinkID", values_to = "value") %>%  # Convert controls to long format
      tidyr::pivot_wider(names_from = Idepic_Bio, values_from = value),  # Transpose controls
    by = "OlinkID"  
  ) %>%
  tibble::column_to_rownames("OlinkID")  

## ICC - https://github.com/courtneyschiffman/Metabolomics-Filtering/blob/master/ICC.R
data_meta_features$ICC <- OmicsProcessing::calculate_ICC(df = data_icc, id_samples = id_samples, id_qc = id_qc)

# distribution plot ====
## Reshape sample data (excluding the first column)
data_features_samples_long <- data_features_samples %>%
  dplyr::select(-1) %>%  # Skip the ID column
  tidyr::pivot_longer(cols = everything(), names_to = "OlinkID", values_to = "value") %>%
  dplyr::mutate(sample_type = "sample")

## Reshape control data (excluding the first column)
data_features_controls_long <- data_features_controls %>%
  dplyr::select(-1) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "OlinkID", values_to = "value") %>%
  dplyr::mutate(sample_type = "control")

## Combine both sample and control data
data_plot <- dplyr::bind_rows(data_features_samples_long, data_features_controls_long)
data_plot <- data_plot %>%
  dplyr::left_join(data_meta_features, by = c("OlinkID"))

## plot distributions with overlaid Shapiro-Wilk W statistics
p <- ggplot2::ggplot(data_plot, 
                     ggplot2::aes(x = value, fill = sample_type)) +
  ggplot2::geom_density(alpha = 0.5) +  # Plot density for each type (Sample/Control)
  ggplot2::facet_wrap(~ OlinkID, 
             scales = "free") +  # Create one plot per protein
  cowplot::theme_cowplot() +
  
  ggplot2::geom_text(ggplot2::aes(label = paste("Sample W:", round(W_stat_samples, 2),
                              "\nControl W:", round(W_stat_controls, 2))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 5) +
  ggplot2::labs(title = "Protein Distributions for Samples and Controls",
       x = "NPX",
       y = "Density",
       fill = "Type") +
  
  ggplot2::scale_fill_manual(values = c("sample" = "blue", "control" = "red")) 

grDevices::tiff(filename = "analysis/001_data-processing/EPIC_olink-immuneonc_distribution-samples.tiff", 
                width = 2000, height = 2000, units = "px")
p  
grDevices::dev.off()

# scatter plot of batch  ====
data_features_samples_long <- data_features_samples %>%
  dplyr::select(Idepic_Bio, everything()) %>%  # Keep the ID column
  tidyr::pivot_longer(cols = -Idepic_Bio, names_to = "OlinkID", values_to = "value") %>%
  dplyr::mutate(sample_type = "sample")

## Join with metadata to get the batch information using Idepic_Bio
data_plot <- data_features_samples_long %>%
  dplyr::left_join(data_meta_samples %>% 
                     dplyr::select(Idepic_Bio, batch_plate, Cncr_Caco_Clrt), 
            by = "Idepic_Bio")  %>%
  dplyr::left_join(data_meta_features, by = "OlinkID")

## order ID by batch and then by cancer status 
id_order <- data_plot %>%
  dplyr::arrange(batch_plate, Cncr_Caco_Clrt, Idepic_Bio) %>%
  dplyr::pull(Idepic_Bio) %>%
  unique()  # Extract unique ordered IDs

data_plot <- data_plot %>%
  dplyr::mutate(Idepic_Bio = factor(Idepic_Bio, levels = id_order))  # Set factor levels based on unique ordered IDs

## plot
p <- ggplot2::ggplot(data_plot, 
                     ggplot2::aes(x = Idepic_Bio, y = value, 
                color = batch_plate,
                shape = as.factor(Cncr_Caco_Clrt))) +
  ggplot2::geom_point(alpha = 0.5) +  # Scatter plot points
  ggplot2::facet_wrap(~ OlinkID, scales = "free") +  # Create one plot per protein
  ggplot2::labs(title = "Scatter Plot of Sample Data by Protein",
       x = "Sample",
       y = "NPX",
       color = "Batch",
       shape = "Case status") +
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()) +
  
  ggplot2::geom_text(ggplot2::aes(label = paste("Missing:", round(missing_pct*100, 2), "%",
                              "\nLOD:", round(LOD, 2),
                              "\nCV:", round(CV, 2),
                              "\nICC:", round(ICC, 2))),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 5, colour = "black") +
  
  ggplot2::geom_hline(ggplot2::aes(yintercept = LOD), linewidth = 2, linetype = "solid", color = "red", )

grDevices::tiff(filename = "analysis/001_data-processing/EPIC_olink-immuneonc_scatter-batch.tiff", 
                width = 3000, height = 3000, units = "px")
p  
grDevices::dev.off()


# write ====
write.table(data_meta_features, "analysis/001_data-processing/EPIC_olink-immuneonc_feature-metadata.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
