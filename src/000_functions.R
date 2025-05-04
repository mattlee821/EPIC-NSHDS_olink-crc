#' Identify Features Below Limit of Detection Threshold
#'
#' This function identifies features where the percentage of values below the
#' limit of detection (LOD) exceeds a specified threshold.
#'
#' @param data_features A dataframe with features as columns and samples as rows.
#' @param data_meta_features A dataframe containing metadata for the features,
#'   including columns for feature identifiers and LOD values.
#' @param feature_col The name of the column in `data_meta_features` that contains
#'   the feature identifiers (which should match the column names in
#'   `data_features`).
#' @param LOD_col The name of the column in `data_meta_features` that contains
#'   the LOD values for each feature.
#' @param percent A numeric value indicating the percentage threshold. Features
#'   with a percentage of values below LOD greater than this value will be
#'   identified.
#'
#' @return A character vector containing the names of the features that exceed the
#'   specified percentage threshold for values below the LOD.
#'
#' @importFrom dplyr select filter %>%
#' @importFrom tibble deframe
#'
#' @examples
#' data_features <- data.frame(
#'   FeatureA = c(10, 5, 2, 8),
#'   FeatureB = c(15, 20, 8, 10)
#' )
#' data_meta <- data.frame(
#'   FeatureID = c("FeatureA", "FeatureB"),
#'   LOD_Value = c(6, 12)
#' )
#'
#' problematic_features <- identify_features_below_lod(
#'   data_features = data_features,
#'   data_meta_features = data_meta,
#'   feature_col = "FeatureID",
#'   LOD_col = "LOD_Value",
#'   percent = 25
#' )
#' print(problematic_features)
#' @export
identify_features_below_lod <- function(
    data_features,
    data_meta_features,
    feature_col,
    LOD_col, 
    percent) {
  
  problematic_features <- character()
  lod_lookup <- data_meta_features %>%
    dplyr::select(dplyr::all_of(feature_col), dplyr::all_of(LOD_col)) %>%
    tibble::deframe() # Convert to a named vector for easy lookup
  
  for (col_name in names(data_features)) {
    if (col_name %in% names(lod_lookup)) {
      lod_value <- lod_lookup[col_name]
      below_lod_count <- data_features %>%
        dplyr::filter(.data[[col_name]] < lod_value) %>%
        nrow()
      percentage_below_lod <- (below_lod_count / nrow(data_features)) * 100
      if (percentage_below_lod > percent) {
        problematic_features <- c(problematic_features, col_name)
      }
    }
  }
  return(problematic_features)
}

#' Calculate Percentage of Values Below Limit of Detection
#'
#' This function calculates the percentage of values below the limit of detection
#' (LOD) for each feature in a given dataset. It handles negative values by
#' treating them as 0 and excludes NA values from the calculation.
#'
#' @param data_features A dataframe with features as columns and samples as rows.
#' @param data_meta_features A dataframe containing metadata for the features,
#'   including columns for feature identifiers and LOD values.
#' @param feature_col The name of the column in `data_meta_features` that contains
#'   the feature identifiers (which should match the column names in
#'   `data_features`). Defaults to "OlinkID".
#' @param LOD_col The name of the column in `data_meta_features` that contains
#'   the LOD values for each feature. Defaults to "LOD".
#'
#' @return A dataframe with two columns:
#'   \itemize{
#'     \item{\code{Feature}:} The name of the feature.
#'     \item{\code{Percentage_Below_LOD}:} The percentage of non-NA values below the LOD for that feature.
#'   }
#'   Features with all NA values will have an `NA` value for
#'   `Percentage_Below_LOD`.
#'
#' @importFrom dplyr select filter %>% mutate bind_rows
#' @importFrom tibble deframe
#'
#' @examples
#' data_features <- data.frame(
#'   FeatureA = c(10, -2, NA, 8),
#'   FeatureB = c(15, 20, 8, NA)
#' )
#' data_meta <- data.frame(
#'   OlinkID = c("FeatureA", "FeatureB"),
#'   LOD = c(6, 12)
#' )
#'
#' lod_percentages <- calculate_below_lod_percentages(
#'   data_features = data_features,
#'   data_meta_features = data_meta
#' )
#' print(lod_percentages)
#' @export
calculate_below_lod_percentages <- function(
    data_features, 
    data_meta_features, 
    feature_col = "OlinkID", 
    LOD_col = "LOD") {

  lod_lookup <- data_meta_features %>%
    dplyr::select(dplyr::all_of(feature_col), dplyr::all_of(LOD_col)) %>%
    tibble::deframe()
  
  results <- data.frame(Feature = character(), Percentage_Below_LOD = numeric())
  
  for (col_name in names(data_features)) {
    if (col_name %in% names(lod_lookup)) {
      lod_value <- lod_lookup[col_name]
      feature_values <- data_features[[col_name]]
      
      # Treat negative values as 0
      feature_values_treated <- ifelse(feature_values < 0, 0, feature_values)
      
      # Exclude NA values
      non_na_values <- feature_values_treated[!is.na(feature_values_treated)]
      
      if (length(non_na_values) > 0) {
        below_lod <- non_na_values < lod_value
        percentage_below_lod <- (sum(below_lod) / length(non_na_values)) * 100
        results <- dplyr::bind_rows(results, data.frame(Feature = col_name, Percentage_Below_LOD = percentage_below_lod))
      } else {
        # Handle cases where all values are NA (optional)
        results <- dplyr::bind_rows(results, data.frame(Feature = col_name, Percentage_Below_LOD = NA_real_))
      }
    }
  }
  return(results)
}

#' Plot Distribution of Features Exceeding Percentage Below LOD Thresholds
#'
#' This function generates a bar chart visualizing the number of features where
#' the calculated percentage of values below the limit of detection (LOD) exceeds
#' specified threshold values.
#'
#' @param df_below_lod_percentages A dataframe, typically the output of
#'   `calculate_below_lod_percentages()`, with columns for Feature and
#'   Percentage_Below_LOD.
#' @param thresholds A numeric vector of threshold values for the
#'   Percentage_Below_LOD. Each threshold will be represented by a bar in the
#'   chart, showing the number of features exceeding that threshold.
#'
#' @return A ggplot object representing the bar chart.
#'
#' @importFrom dplyr filter %>%
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme element_text
#' @importFrom cowplot theme_cowplot
#'
#' @examples
#' data_features <- data.frame(
#'   FeatureA = c(10, -2, NA, 8),
#'   FeatureB = c(15, 20, 8, NA),
#'   FeatureC = c(0.5, 1.2, -0.3, NA)
#' )
#' data_meta <- data.frame(
#'   OlinkID = c("FeatureA", "FeatureB", "FeatureC"),
#'   LOD = c(6, 12, 1)
#' )
#'
#' lod_percentages <- calculate_below_lod_percentages(
#'   data_features = data_features,
#'   data_meta_features = data_meta
#' )
#'
#' plot_features_exceeding_threshold(
#'   df_below_lod_percentages = lod_percentages,
#'   thresholds = c(10, 50)
#' )
#' @export
plot_features_exceeding_threshold <- function(
    df_below_lod_percentages, 
    thresholds = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100)) {
  
  results_exceeding <- data.frame(Threshold = numeric(), Number_of_Features = integer())
  
  for (threshold in thresholds) {
    num_exceeding <- df_below_lod_percentages %>%
      dplyr::filter(Percentage_Below_LOD > threshold) %>%
      nrow()
    results_exceeding <- dplyr::bind_rows(results_exceeding, data.frame(Threshold = threshold, Number_of_Features = num_exceeding))
  }
  
  ggplot2::ggplot(results_exceeding, 
                  ggplot2::aes(x = factor(Threshold), 
                               y = Number_of_Features)) +
    ggplot2::geom_bar(stat = "identity", fill = "purple") +
    ggplot2::geom_text(ggplot2::aes(label = Number_of_Features), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "N features with > X% of\n samples with values < LOD",
      x = "Threshold percentage",
      y = "N features"
    ) +
    cowplot::theme_cowplot() 
}

#' Identify Samples Below Limit of Detection Threshold
#'
#' This function identifies samples where the percentage of features below the
#' limit of detection (LOD) exceeds a specified threshold.
#'
#' @param data_features A dataframe with features as columns and samples as rows.
#'   The column specified by `sample_col` should contain the sample IDs.
#' @param data_meta_features A dataframe containing metadata for the features,
#'   including columns for feature identifiers and LOD values.
#' @param feature_col The name of the column in `data_meta_features` that
#'   contains the feature identifiers (which should match the column names in
#'   `data_features`).
#' @param LOD_col The name of the column in `data_meta_features` that contains
#'   the LOD values for each feature.
#' @param sample_col The name of the column in `data_features` that contains
#'   the sample IDs.
#' @param percent A numeric value indicating the percentage threshold. Samples
#'   with a percentage of features below LOD greater than this value will be
#'   identified.
#'
#' @return A character vector containing the sample IDs that exceed the
#'   specified percentage threshold for features below the LOD.
#'
#' @importFrom dplyr select filter %>%
#' @importFrom tibble deframe column_to_rownames
#'
#' @examples
#' data_features <- data.frame(
#'   SampleID = c("Sample1", "Sample2", "Sample3", "Sample4"),
#'   FeatureA = c(10, 5, 2, 8),
#'   FeatureB = c(15, 20, 8, 10)
#' )
#' data_meta <- data.frame(
#'   FeatureID = c("FeatureA", "FeatureB"),
#'   LOD_Value = c(6, 12)
#' )
#'
#' problematic_samples <- identify_samples_below_lod(
#'   data_features = data_features,
#'   data_meta_features = data_meta,
#'   feature_col = "FeatureID",
#'   LOD_col = "LOD_Value",
#'   sample_col = "SampleID",
#'   percent = 25
#' )
#' print(problematic_samples)
#' @export
identify_samples_below_lod <- function(
    data_features,
    data_meta_features,
    feature_col,
    LOD_col,
    sample_col,
    percent) {
  
  problematic_samples <- character()
  lod_lookup <- data_meta_features %>%
    dplyr::select(dplyr::all_of(feature_col), dplyr::all_of(LOD_col)) %>%
    tibble::deframe()
  
  # Convert the sample ID column to row names
  data_features_with_rownames <- data_features %>%
    tibble::column_to_rownames(var = sample_col)
  
  for (sample_name in rownames(data_features_with_rownames)) {
    sample_values <- data_features_with_rownames[sample_name, ]
    below_lod_count <- 0
    total_feature_count <- 0
    
    for (feature in names(sample_values)) {
      if (feature %in% names(lod_lookup)) {
        lod_value <- lod_lookup[feature]
        if (!is.na(sample_values[[feature]])) {
          total_feature_count <- total_feature_count + 1
          if (sample_values[[feature]] < lod_value) {
            below_lod_count <- below_lod_count + 1
          }
        }
      }
    }
    
    # Avoid division by zero if no features have non-NA values
    if (total_feature_count > 0) {
      percentage_below_lod <- (below_lod_count / total_feature_count) * 100
    } else {
      percentage_below_lod <- 0  # Or NA, depending on desired behavior
    }
    
    if (percentage_below_lod > percent) {
      problematic_samples <- c(problematic_samples, sample_name)
    }
  }
  return(problematic_samples)
}

#' Calculate Percentage of Features Below Limit of Detection for Samples
#'
#' This function calculates the percentage of features below the limit of
#' detection (LOD) for each sample in a given dataset. It handles negative
#' values by treating them as 0 and excludes NA values from the calculation.
#'
#' @param data_features A dataframe with features as columns and samples as rows.
#'   The first column should contain the sample IDs.
#' @param data_meta_features A dataframe containing metadata for the features,
#'   including columns for feature identifiers and LOD values.
#' @param feature_col The name of the column in `data_meta_features` that
#'   contains the feature identifiers (which should match the column names in
#'   `data_features`). Defaults to "OlinkID".
#' @param LOD_col The name of the column in `data_meta_features` that contains
#'   the LOD values for each feature. Defaults to "LOD".
#' @param sample_col The name of the column in `data_features` that contains
#'   the sample IDs. Defaults to the first column.
#'
#' @return A dataframe with two columns:
#'   \itemize{
#'     \item{\code{Sample}:} The sample IDs.
#'     \item{\code{Percentage_Below_LOD}:} The percentage of non-NA features
#'       below the LOD for that sample.
#'   }
#'   Samples with all NA values will have an `NA` value for
#'   `Percentage_Below_LOD`.
#'
#' @importFrom dplyr select filter %>% mutate bind_rows
#' @importFrom tibble deframe
#'
#' @examples
#' data_features <- data.frame(
#'   SampleID = c("Sample1", "Sample2", "Sample3", "Sample4"),
#'   FeatureA = c(10, -2, NA, 8),
#'   FeatureB = c(15, 20, 8, NA)
#' )
#' data_meta <- data.frame(
#'   OlinkID = c("FeatureA", "FeatureB"),
#'   LOD = c(6, 12)
#' )
#'
#' lod_percentages <- calculate_below_lod_percentages_samples(
#'   data_features = data_features,
#'   data_meta_features = data_meta,
#'   feature_col = "OlinkID",
#'   LOD_col = "LOD",
#'   sample_col = "SampleID"
#' )
#' print(lod_percentages)
#' @export
calculate_below_lod_percentages_samples <- function(
    data_features,
    data_meta_features,
    feature_col = "OlinkID",
    LOD_col = "LOD",
    sample_col = "Idepic_Bio") {
  
  lod_lookup <- data_meta_features %>%
    dplyr::select(dplyr::all_of(feature_col), dplyr::all_of(LOD_col)) %>%
    tibble::deframe()
  
  # Convert the sample ID column to row names
  data_features_with_rownames <- data_features %>%
    tibble::column_to_rownames(var = sample_col)
  
  results <- data.frame(Sample = character(), Percentage_Below_LOD = numeric())
  
  for (sample_name in rownames(data_features_with_rownames)) {
    sample_values <- data_features_with_rownames[sample_name, ]
    feature_count <- 0
    below_lod_count <- 0
    
    for (feature in names(sample_values)) {
      if (feature %in% names(lod_lookup)) {
        lod_value <- lod_lookup[feature]
        if (!is.na(sample_values[[feature]])) {
          feature_count <- feature_count + 1
          # Treat negative values as 0
          value_treated <- ifelse(sample_values[[feature]] < 0, 0, sample_values[[feature]])
          if (value_treated < lod_value) {
            below_lod_count <- below_lod_count + 1
          }
        }
      }
    }
    
    if (ncol(data_features_with_rownames) > 0) { # Check if there are features
      percentage_below_lod <- (below_lod_count / ncol(data_features_with_rownames)) * 100
    } else {
      percentage_below_lod <- NA_real_ # If there are no features
    }
    results <- dplyr::bind_rows(results, data.frame(Sample = sample_name, Percentage_Below_LOD = percentage_below_lod))
  }
  return(results)
}

#' Plot Distribution of Samples Exceeding Percentage Below LOD Thresholds
#'
#' This function generates a bar chart visualizing the number of samples where
#' the calculated percentage of features below the limit of detection (LOD)
#' exceeds specified threshold values.
#'
#' @param df_below_lod_percentages A dataframe, typically the output of
#'   `calculate_below_lod_percentages_samples()`, with columns for Sample and
#'   Percentage_Below_LOD.
#' @param thresholds A numeric vector of threshold values for the
#'   Percentage_Below_LOD. Each threshold will be represented by a bar in the
#'   chart, showing the number of samples exceeding that threshold.
#' @param sample_col The name of the column in `df_below_lod_percentages` that
#'   contains the sample IDs. Defaults to "Sample".
#'
#' @return A ggplot object representing the bar chart.
#'
#' @importFrom dplyr filter %>%
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme element_text
#' @importFrom cowplot theme_cowplot
#'
#' @examples
#' data_features <- data.frame(
#'   SampleID = c("Sample1", "Sample2", "Sample3", "Sample4"),
#'   FeatureA = c(10, -2, NA, 8),
#'   FeatureB = c(15, 20, 8, NA),
#'   FeatureC = c(0.5, 1.2, -0.3, NA)
#' )
#' data_meta <- data.frame(
#'   OlinkID = c("FeatureA", "FeatureB", "FeatureC"),
#'   LOD = c(6, 12, 1)
#' )
#'
#' lod_percentages <- calculate_below_lod_percentages_samples(
#'   data_features = data_features,
#'   data_meta_features = data_meta,
#'   feature_col = "OlinkID",
#'   LOD_col = "LOD",
#'   sample_col = "SampleID"
#' )
#'
#' plot_samples_exceeding_threshold <- plot_samples_exceeding_threshold(
#'   df_below_lod_percentages = lod_percentages,
#'   thresholds = c(10, 50),
#'   sample_col = "Sample"
#' )
#' @export
plot_samples_exceeding_threshold <- function(
    df_below_lod_percentages,
    thresholds = c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100),
    sample_col = "Sample") {
  
  results_exceeding <- data.frame(Threshold = numeric(), Number_of_Samples = integer())
  
  for (threshold in thresholds) {
    num_exceeding <- df_below_lod_percentages %>%
      dplyr::filter(Percentage_Below_LOD > threshold) %>%
      nrow()
    results_exceeding <- dplyr::bind_rows(results_exceeding, data.frame(Threshold = threshold, Number_of_Samples = num_exceeding))
  }
  
  ggplot2::ggplot(results_exceeding,
                  ggplot2::aes(x = factor(Threshold), y = Number_of_Samples)) +
    ggplot2::geom_bar(stat = "identity", fill = "purple") +
    ggplot2::geom_text(ggplot2::aes(label = Number_of_Samples), vjust = -0.5, size = 3) +
    ggplot2::labs(
      title = "N samples with > X% of\n samples with values < LOD",
      x = "Threshold percentage",
      y = "N samples"
    ) +
    cowplot::theme_cowplot() 
}

#' Analyze Data with Multiple Models
#'
#' This function iterates through a list of datasets, applies a set of statistical models,
#' and extracts and returns the results in a data frame.
#'
#' @param data_list A list of data frames, where each data frame represents a dataset to be analyzed.
#' @param models A list of functions, where each function represents a statistical model to be applied.
#' @param covariates (Optional) A character string specifying the names of covariates to be used in 'model_2'.
#'                     If NULL, a default set of covariates is used.
#'                     Defaults to NULL.
#'
#' @return A data frame containing the extracted results from each model applied to each dataset.
#'         The columns include:
#'         \itemize{
#'           \item \code{exposure}: The name of the exposure variable.
#'           \item \code{outcome}: The outcome variable (extracted from the dataset name).
#'           \item \code{sex}: Sex information (extracted from the dataset name).
#'           \item \code{model}: The name of the model applied.
#'           \item \code{followup}: Follow-up information (extracted from the dataset name).
#'           \item \code{n}: Number of observations.
#'           \item \code{nevent}: Number of events (if applicable).
#'           \item \code{coef}: Coefficient estimate for the exposure variable.
#'           \item \code{coef_exp}: Exponentiated coefficient estimate (e.g., odds ratio, hazard ratio).
#'           \item \code{se}: Standard error of the coefficient.
#'           \item \code{ci_lower_exp}: Lower confidence interval for the exponentiated coefficient.
#'           \item \code{ci_upper_exp}: Upper confidence interval for the exponentiated coefficient.
#'           \item \code{se_robust}: Robust standard error of the coefficient (if available).
#'           \item \code{pval}: P-value for the coefficient.
#'           \item \code{degrees_freedom}: Degrees of freedom for the model test.
#'           \item \code{concordance}: Concordance statistic (if applicable).
#'           \item \code{concordance_se}: Standard error of the concordance statistic.
#'           \item \code{likelihood_ratio_test}: Likelihood ratio test statistic.
#'           \item \code{likelihood_ratio_test_pval}: P-value for the likelihood ratio test.
#'           \item \code{wald_test}: Wald test statistic.
#'           \item \code{wald_test_pval}: P-value for the Wald test.
#'           \item \code{score_test}: Score test statistic.
#'           \item \code{score_test_pval}: P-value for the Score test.
#'           \item \code{robust_score_test_pval}: P-value for the robust score test.
#'           \item \code{exclusion_missing}: Number of observations excluded due to missingness.
#'         }
#'
#' @details
#' This function automates the process of applying multiple statistical models to a
#' series of datasets. It iterates through each dataset in the input list, identifies
#' potential exposure variables, and then applies each model in the 'models' list.
#' The results from each model are extracted and combined into a single data frame.
#'
#' The function handles potential errors during model fitting and extracts key
#' statistics from the fitted model objects. It also provides a mechanism to specify
#' covariates for a particular model ('model_2').
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#'
#' @examples
#' # Example usage (assuming 'my_data_list' and 'my_models' are defined)
#' # results_df <- analysis(data_list = my_data_list, models = my_models)
#' # results_df_with_covariates <- analysis(data_list = my_data_list, models = my_models, covariates = "age + sex")
#' # print(results_df)
#' # print(results_df_with_covariates)
#'
#' @export
analysis <- function(data_list, models, covariates = NULL) {
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each dataset in the input list
  for (i in seq_along(data_list)) {
    dataset_name <- names(data_list[i])
    dataset <- data_list[[i]]
    
    print(paste("# Analyzing:", dataset_name))
    
    # Iterate through potential exposure variables
    exposure_vars <- names(dataset)[grepl("OID", names(dataset))]
    
    for (exposure_var in exposure_vars) {
      # Check for missing values in exposure variable
      if (any(!is.na(dataset[[exposure_var]]))) {
        # Iterate through models
        for (model_name in names(models)) {
          cat(paste0("      + ", model_name, "; exposure: ", exposure_var, "; data: ", dataset_name, "\n"))
          
          # Call the appropriate model function - if covariates are not provided stop 
          if (model_name == "model_2") {
            if (is.null(covariates) || length(covariates) == 0 || all(is.na(covariates)) || all(covariates == "")) {
              stop("model_2 requires covariates; please provide covariates")
            }
            # covariates check passed, assign the covariates string
            covariates_string <- covariates
          }
          
          tryCatch({
            if (model_name == "model_2") {
              model_result <- models[[model_name]](dataset, exposure_var, covariates = covariates_string)
            } else {
              model_result <- models[[model_name]](dataset, exposure_var)
            }
            
            # Extract model results (using helper function)
            result_row <- extract_model_results(model_result)
            result_row$exposure <- exposure_var
            result_row$outcome <- strsplit(dataset_name, "_")[[1]][1]
            result_row$sex <- strsplit(dataset_name, "_")[[1]][2]
            result_row$model <- model_name
            result_row$followup <- strsplit(dataset_name, "_")[[1]][3]
            
            results_list <- append(results_list, list(result_row))
            
          }, error = function(e) {
            # Handle potential errors gracefully
            warning(paste("Error in model", model_name, "for exposure", exposure_var, ":", e$message))
            # Add a row of NA values in case of error
            results_list <- append(results_list, list(data.frame(
              exposure = exposure_var,
              outcome = strsplit(dataset_name, "_")[[1]][1],
              sex = strsplit(dataset_name, "_")[[1]][2],
              model = model_name,
              followup = strsplit(dataset_name, "_")[[1]][3],
              n = NA,
              nevent = NA,
              exclusion_missing = NA,
              coef = NA,
              coef_exp = NA,
              se = NA,
              ci_lower_exp = NA,
              ci_upper_exp = NA,
              ci_upper_exp = NA,
              se_robust = NA,
              pval = NA,
              degrees_freedom = NA,
              concordance = NA,
              concordance_se = NA,
              likelihood_ratio_test = NA,
              likelihood_ratio_test_pval = NA,
              wald_test = NA,
              wald_test_pval = NA,
              score_test = NA,
              score_test_pval = NA
              # robust_score_test_pval = NA,
            )))
          })
        }
      } else {
        # Add a row of NA values if exposure is completely missing
        results_list <- append(results_list, list(data.frame(
          exposure = exposure_var,
          outcome = strsplit(dataset_name, "_")[[1]][1],
          sex = strsplit(dataset_name, "_")[[1]][2],
          model = NA,
          followup = strsplit(dataset_name, "_")[[1]][3],
          n = NA,
          nevent = NA,
          exclusion_missing = NA,
          coef = NA,
          coef_exp = NA,
          se = NA,
          ci_lower_exp = NA,
          ci_upper_exp = NA,
          ci_upper_exp = NA,
          se_robust = NA,
          pval = NA,
          degrees_freedom = NA,
          concordance = NA,
          concordance_se = NA,
          likelihood_ratio_test = NA,
          likelihood_ratio_test_pval = NA,
          wald_test = NA,
          wald_test_pval = NA,
          score_test = NA,
          # robust_score_test_pval = NA
        )))
      }
    }
  }
  
  # Convert the list of results to a dataframe
  data_result <- do.call(rbind, results_list)
  data_result <- data_result %>%
    dplyr::select(
      exposure, outcome, sex, model, followup,
      everything()
    )
  return(data_result)
}

#' Extract Model Results
#'
#' This function extracts key results from a fitted statistical model object.
#'
#' @param model_result A fitted model object (e.g., from `clogit`, `coxph`, `glm`).
#' @return A data frame with one row containing extracted model results.
#'         The columns include:
#'         \itemize{
#'           \item \code{n}: Number of observations.
#'           \item \code{nevent}: Number of events (if applicable).
#'           \item \code{coef}: Coefficient estimate for the exposure variable.
#'           \item \code{coef_exp}: Exponentiated coefficient estimate (e.g., odds ratio, hazard ratio).
#'           \item \code{se}: Standard error of the coefficient.
#'           \item \code{ci_lower_exp}: Lower confidence interval for the exponentiated coefficient.
#'           \item \code{ci_upper_exp}: Upper confidence interval for the exponentiated coefficient.
#'           \item \code{pval}: P-value for the coefficient.
#'           \item \code{degrees_freedom}: Degrees of freedom for the model test.
#'           \item \code{concordance}: Concordance statistic (if applicable).
#'           \item \code{concordance_se}: Standard error of the concordance statistic.
#'           \item \code{likelihood_ratio_test}: Likelihood ratio test statistic.
#'           \item \code{likelihood_ratio_test_pval}: P-value for the likelihood ratio test.
#'           \item \code{wald_test}: Wald test statistic.
#'           \item \code{wald_test_pval}: P-value for the Wald test.
#'           \item \code{score_test}: Score test statistic.
#'           \item \code{score_test_pval}: P-value for the Score test.
#'           \item \code{exclusion_missing}: Number of observations excluded due to missingness.
#'         }
#'
#' @details
#' This function attempts to extract relevant statistics from a fitted model object.
#' It uses `tryCatch` to handle potential errors during the extraction process.
#' The extracted statistics may vary depending on the type of model fitted.
#'
extract_model_results <- function(model_result) {
  # Initialize a data frame with NA values
  result_row <- data.frame(
    n = NA,
    nevent = NA,
    exclusion_missing = NA,
    coef = NA,
    coef_exp = NA,
    se = NA,
    ci_lower_exp = NA,
    ci_upper_exp = NA,
    # se_robust = NA,
    pval = NA,
    degrees_freedom = NA,
    concordance = NA,
    concordance_se = NA,
    likelihood_ratio_test = NA,
    likelihood_ratio_test_pval = NA,
    wald_test = NA,
    wald_test_pval = NA,
    score_test = NA,
    score_test_pval = NA
    # robust_score_test_pval = NA
  )
  
  tryCatch({
    # Extract values from model_result
    s <- summary(model_result)
    
    if (!is.null(s$coef) && nrow(s$coef) > 0) {
      result_row$n <- ifelse(is.numeric(s$n), s$n, NA)
      result_row$nevent <- ifelse(is.numeric(s$nevent), s$nevent, NA)
    }
    
    if (!is.null(model_result$na.action)) {
      result_row$exclusion_missing <- ifelse(is.numeric(length(model_result$na.action)), length(model_result$na.action), NA)
    }
    
    if (!is.null(s$coef) && nrow(s$coef) > 0) {
      result_row$coef <- ifelse(is.numeric(s$coef[1, 1]), s$coef[1, 1], NA)
      result_row$se <- ifelse(is.numeric(s$coef[1, "se(coef)"]), s$coef[1, "se(coef)"], NA)
      # result_row$se_robust <- ifelse(is.numeric(s$coef[1, "robust se"]), s$coef[1, "robust se"], NA)
      result_row$pval <- ifelse(is.numeric(s$coef[1, "Pr(>|z|)"]), s$coef[1, "Pr(>|z|)"], NA)
    }
    
    if (!is.null(s$conf.int) && nrow(s$conf.int) > 0) {
      result_row$coef_exp <- ifelse(is.numeric(s$conf.int[1, 1]), s$conf.int[1, 1], NA)
      result_row$ci_lower_exp <- ifelse(is.numeric(s$conf.int[1, 3]), s$conf.int[1, 3], NA)
      result_row$ci_upper_exp <- ifelse(is.numeric(s$conf.int[1, 4]), s$conf.int[1, 4], NA)
    }
    
    if (!is.null(s$logtest)) {
      result_row$degrees_freedom <- ifelse(is.numeric(s$logtest[["df"]]), s$logtest[["df"]], NA)
      result_row$likelihood_ratio_test <- ifelse(is.numeric(s$logtest[["test"]]), s$logtest[["test"]], NA)
      result_row$likelihood_ratio_test_pval <- ifelse(is.numeric(s$logtest[["pvalue"]]), s$logtest[["pvalue"]], NA)
    }
    
    if (!is.null(s$waldtest)) {
      result_row$wald_test <- ifelse(is.numeric(s$waldtest[["test"]]), s$waldtest[["test"]], NA)
      result_row$wald_test_pval <- ifelse(is.numeric(s$waldtest[["pvalue"]]), s$waldtest[["pvalue"]], NA)
    }
    
    if (!is.null(s$sctest)) {
      result_row$score_test <- ifelse(is.numeric(s$sctest[["test"]]), s$sctest[["test"]], NA)
      result_row$score_test_pval <- ifelse(is.numeric(s$sctest[["pvalue"]]), s$sctest[["pvalue"]], NA)
    }
    
    # if (!is.null(s$robscore)) {
    #   result_row$robust_score_test <- ifelse(is.numeric(s$robscore[["test"]]), s$robscore[["test"]], NA)
    #   result_row$robust_score_test_pval <- ifelse(is.numeric(s$robscore[["pvalue"]]), s$robscore[["pvalue"]], NA)
    # }
    
    if (!is.null(model_result$concordance)) {
      result_row$concordance <- ifelse(is.numeric(model_result$concordance[["concordance"]]), model_result$concordance[["concordance"]], NA)
      result_row$concordance_se <- ifelse(is.numeric(model_result$concordance[["std"]]), model_result$concordance[["std"]], NA)
    }
    
  }, error = function(e) {
    warning(paste("Error extracting results:", e$message))
  })
  
  return(result_row)
}

#' List of Statistical Models
#'
#' This list contains functions for fitting different statistical models.
#'
#' @format A list of two functions:
#'   \describe{
#'     \item{\code{model_1}}{A function to fit a conditional logistic regression model with a specific set of covariates.}
#'     \item{\code{model_2}}{A function to fit a conditional logistic regression model with user-specified covariates.}
#'   }
#' @details
#'   The functions in this list are designed to be used in epidemiological analyses,
#'   particularly for case-control studies. They utilize the `survival::clogit`
#'   function for conditional logistic regression.
#'
#' @name models
#' @rdname models
#' @export
#'
#' @examples
#' # Example usage (assuming 'df' and 'exposure_var' are defined)
#' # results_model_1 <- models$model_1(df = df, exposure_var = "my_exposure")
#' # results_model_2 <- models$model_2(df = df, exposure_var = "my_exposure", covariates = "age + sex")
#' # print(results_model_1)
#' # print(results_model_2)
#'
models <- list(
  #' Conditional Logistic Regression Model 1
  #'
  #' Fits a conditional logistic regression model with a pre-defined set of covariates.
  #'
  #' @param df A data frame containing the data for the analysis.
  #' @param exposure_var A character string specifying the name of the exposure variable in `df`.
  #' @return A fitted `clogit` model object (of class `coxph`).
  #' @import survival
  #' @export
  model_1 = function(df, exposure_var){
    survival::clogit(Cncr_Caco_Clrt ~ get(exposure_var) +
                       Fasting_C +
                       day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi +
                       survival::strata(Match_Caseset),
                     data = df,
                     method = 'exact',
                     na.action = 'na.exclude')
  },
  
  #' Conditional Logistic Regression Model 2
  #'
  #' Fits a conditional logistic regression model with user-specified covariates.
  #'
  #' @param df A data frame containing the data for the analysis.
  #' @param exposure_var A character string specifying the name of the exposure variable in `df`.
  #' @param covariates A character string specifying the names of the covariates to include in the model,
  #'                   separated by '+'. For example: "age + sex + bmi".
  #' @return A fitted `clogit` model object (of class `coxph`).
  #' @import survival
  #' @export
  model_2 = function(df, exposure_var, covariates){
    formula_str <- paste0(
      "Cncr_Caco_Clrt ~ ", exposure_var, " + ",
      "Fasting_C + ",
      "day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_cos_4_pi + day_of_year_blood_draw_sin_2_pi + day_of_year_blood_draw_sin_4_pi + ",
      covariates, " + ",
      "survival::strata(Match_Caseset)"
    )
    
    survival::clogit(as.formula(formula_str),
                     data = df,
                     method = 'exact',
                     na.action = 'na.exclude')
  }
)

