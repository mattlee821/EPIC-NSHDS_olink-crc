rm(list=ls())
set.seed(821)

# environment ====
source("src/000_source.R")
rm(list=ls())

# data_samples ====
data_samples <- haven::read_sas(data_file = "data/raw/EPIC_phenotype/clrt_caco.sas7bdat")
data_samples <- data_samples %>%
  # make age groups
  dplyr::mutate(
    age_group = as.integer(cut(
      Age_Recr,
      breaks = seq(0, max(Age_Recr, na.rm = TRUE) + 5, by = 5),
      right = FALSE,
      labels = seq(0, max(Age_Recr, na.rm = TRUE), by = 5)
    )),
    
    # make day of blood draw
    day_of_year = data.table::yday(as.Date(D_Bld_Coll, format = "%Y-%m-%d")), # Convert date of blood collection to day of year
    ## Trigonometric transformations
    ## Regarding adjustment for month/season of blood draw, in analyses of circulating vitamin D we have had success with spherical adjustment for day-of-year of blood draw. 
    ## By using a (very) truncated Fourier transformation, I think we will be able to get away with four parameters for this adjustment rather than the 11 parameters needed to adjust for month. 
    ## This also has the added bonus that the adjustment is smooth, with no artificial discontinuities from season to season or year to year.
    day_of_year_blood_draw_sin_2_pi = sin(2 * pi * day_of_year / 365),
    day_of_year_blood_draw_cos_2_pi = cos(2 * pi * day_of_year / 365),
    day_of_year_blood_draw_sin_4_pi = sin(4 * pi * day_of_year / 365),
    day_of_year_blood_draw_cos_4_pi = cos(4 * pi * day_of_year / 365),
    
    # make fasting status
    fasting_c_misscat = dplyr::case_when(
      Fasting_C == "No" ~ "0", 
      Fasting_C == "In between" ~ "1",
      Fasting_C == "Yes" ~ "2",
      is.na(Fasting_C) ~ "Missing"
    ),
    
    # covariate formatting
    # Smoking status: 1 = never; 2 = former; 3 = smoker; 4 = unknown
    Smoke_Stat = dplyr::case_when(
      is.na(Smoke_Stat) ~ "Unknown",
      Smoke_Stat == 1 ~ "Never",
      Smoke_Stat == 2 ~ "Former",
      Smoke_Stat == 3 ~ "Smoker",
      Smoke_Stat == 4 ~ "Unknown"
    ) %>% factor(levels = c("Never", "Former", "Smoker", "Unknown")),
    # Physical activity index: 1 = inactive; 2 = moderately inactive; 3 = moderately active; 4 = active; 5 = missing
    Pa_Index = dplyr::case_when(
      is.na(Pa_Index) ~ "Missing",
      Pa_Index == 1 ~ "Inactive",
      Pa_Index == 2 ~ "Moderately inactive",
      Pa_Index == 3 ~ "Moderately active",
      Pa_Index == 4 ~ "Active",
      Pa_Index == 5 ~ "Missing"
    ) %>% factor(levels = c("Inactive", "Moderately inactive", "Moderately active", "Active", "Missing")),
    # Level of schooling: 0 = none; 1 = primary; 2 = technical/professional; 3 = secondary; 4 = longer; 5 = not specified
    L_School = dplyr::case_when(
      is.na(L_School) ~ "Not specified",
      L_School == 0 ~ "None",
      L_School == 1 ~ "Primary",
      L_School == 2 ~ "Technical/professional",
      L_School == 3 ~ "Secondary",
      L_School == 4 ~ "Longer education",
      L_School == 5 ~ "Not specified"
    ) %>% factor(levels = c("None", "Primary", "Secondary", "Technical/professional", "Longer education", "Not specified", "Missing"))
  )

temp_mapping <- data_samples %>%
  dplyr::select(Idepic_Bio) %>%
  tidyr::extract(Idepic_Bio, into = c("prefix", "SampleID"), regex = "^([^_]+)_(.+)$", remove = FALSE) %>%
  dplyr::mutate(SampleID = gsub("_", "", SampleID)) %>%
  dplyr::select(-prefix)

# data olink ====
data_immuneonc <- OlinkAnalyze::read_NPX(filename = "data/raw/EPIC_olink-immuneonc/20211440_Murphy_NPX.xlsx") %>%
  tibble::as_tibble()
data_explore <- OlinkAnalyze::read_NPX(filename = "data/raw/EPIC_olink-explore/WL-3904_NPX_2024-12-09.csv") %>%
  tibble::as_tibble()

# formatting ====
data_immuneonc <- data_immuneonc %>%
  dplyr::mutate(Sample_Type = ifelse(grepl("CONTROL", SampleID), "CONTROL", "SAMPLE"))
data_explore <- data_explore %>%
  dplyr::mutate(Sample_Type = ifelse(SampleID == "Empty well", "CONTROL", Sample_Type)) %>%
  dplyr::mutate(
    SampleID = if_else(
      Sample_Type == "SAMPLE",
      stringr::str_pad(as.character(as.numeric(SampleID)), width = 5, pad = "0"),
      SampleID
    )
  ) %>%
  dplyr::left_join(temp_mapping, by = "SampleID") %>%
  dplyr::mutate(SampleID = ifelse(Sample_Type == "SAMPLE", Idepic_Bio, SampleID)) %>%
  dplyr::select(-Idepic_Bio)

# map Assay to explore
data_info_immuneonc <- data_immuneonc %>%
  dplyr::select(OlinkID, Assay, UniProt) %>%
  unique()
data_info_explore <- data_explore %>%
  dplyr::select(OlinkID, Assay, UniProt) %>%
  unique()

data_info <- data_info_explore %>%
  dplyr::left_join(data_info_immuneonc, by = "UniProt") %>%
  dplyr::rename(Assay_explore = Assay.x,
                Assay_immuneonc = Assay.y,
                OlinkID_explore = OlinkID.x,
                OlinkID_immuneonc = OlinkID.y) %>%
  # Create unified OlinkID column for indexing
  dplyr::mutate(OlinkID = coalesce(OlinkID_explore, OlinkID_immuneonc),
                Assay_base = coalesce(Assay_explore, Assay_immuneonc)) %>%
  dplyr::group_by(Assay_base) %>%
  dplyr::mutate(index = row_number(),
                Assay_map = if (n() > 1) paste0(Assay_base, ";", index) else Assay_base) %>%
  dplyr::ungroup() %>%
  dplyr::select(-OlinkID, -Assay_base, -index)

data_info_immuneonc <- data_info %>%
  dplyr::select(OlinkID_immuneonc, UniProt, Assay_explore) %>%
  tidyr::drop_na() %>%
  unique() %>%
  dplyr::rename(Assay = Assay_explore)

data_info_explore <- data_info %>%
  dplyr::select(OlinkID_explore, UniProt, Assay_map) %>%
  tidyr::drop_na() %>%
  unique() %>%
  dplyr::rename(Assay = Assay_map)

data_immuneonc <- data_immuneonc %>%
  dplyr::left_join(data_info_immuneonc, by = c("OlinkID" = "OlinkID_immuneonc", "UniProt")) %>%
  dplyr::mutate(Assay = Assay.y) %>%  
  dplyr::select(-Assay.x, -Assay.y)   

data_explore <- data_explore %>%
  dplyr::left_join(data_info_explore, by = c("OlinkID" = "OlinkID_explore", "UniProt")) %>%
  dplyr::mutate(Assay = Assay.y) %>%
  dplyr::select(-Assay.x, -Assay.y) 

# pheno join ====
data_immuneonc <- data_immuneonc %>%
  dplyr::left_join(data_samples, by = c("SampleID" = "Idepic_Bio"))
data_explore <- data_explore %>%
  dplyr::left_join(data_samples, by = c("SampleID" = "Idepic_Bio"))

# exclusions ====
## exclusion 1 immuneonc: exclude samples present in multiple casesets but which are on different plates/batches ====
# Step 1: Identify problematic Match_Casesets with > 2 unique SampleIDs
problematic_casesets <- data_immuneonc %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::summarise(n_SampleIDs = n_distinct(SampleID), .groups = "drop") %>%
  dplyr::filter(n_SampleIDs > 2) %>%
  dplyr::pull(Match_Caseset)
# Step 2: Filter the main data to those problematic casesets,
# keeping relevant columns: SampleID, Match_Caseset, PlateID, Cncr_Caco_Clrt
problematic_data <- data_immuneonc %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::filter(Match_Caseset %in% problematic_casesets) %>%
  dplyr::select(SampleID, Match_Caseset, PlateID, Cncr_Caco_Clrt) %>%
  dplyr::distinct()
# Step 3: Identify Match_Casesets that have multiple PlateIDs (potential problem)
casesets_multiple_plates <- problematic_data %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::summarise(n_PlateIDs = n_distinct(PlateID), .groups = "drop") %>%
  dplyr::filter(n_PlateIDs > 1) %>%
  dplyr::pull(Match_Caseset)
# Step 4: Filter data to only these problematic Match_Casesets with multiple PlateIDs
problematic_multi_plate <- problematic_data %>%
  dplyr::filter(Match_Caseset %in% casesets_multiple_plates)
# Step 5: Find SampleIDs that are the *only* sample on their PlateID within that Match_Caseset (candidates for removal)
removal_candidates <- problematic_multi_plate %>%
  dplyr::group_by(Match_Caseset, PlateID) %>%
  dplyr::filter(n() == 1) %>%
  dplyr::ungroup()
# Step 6: Define a function to check if removing a sample keeps both cases and controls in the Match_Caseset
check_keep_case_control <- function(caseset, sample_to_remove, data) {
  remaining <- data %>%
    dplyr::filter(Match_Caseset == caseset, SampleID != sample_to_remove)
  
  has_case <- any(remaining$Cncr_Caco_Clrt == 1)
  has_control <- any(remaining$Cncr_Caco_Clrt == 0)
  
  return(has_case & has_control)
}
# Step 7: Apply the function to each candidate to test if they can be safely removed
removal_candidates <- removal_candidates %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    can_remove = check_keep_case_control(Match_Caseset, SampleID, problematic_data)
  ) %>%
  dplyr::ungroup()
# Step 8: Filter to samples that can actually be removed safely
removable_samples <- removal_candidates %>%
  dplyr::filter(can_remove)
# Step 9: Output removable samples, or remove them from original data if desired
VAR_sample_exclude <- removable_samples %>%
  dplyr::select(SampleID, Match_Caseset, PlateID)
# Step 10: exclude
data_immuneonc <- data_immuneonc %>%
  dplyr::anti_join(VAR_sample_exclude, by = c("SampleID", "Match_Caseset", "PlateID"))

## exclusion 1 explore: exclude samples present in multiple casesets but which are on different plates/batches ====
# Step 1: Identify problematic Match_Casesets with > 2 unique SampleIDs
problematic_casesets <- data_explore %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::summarise(n_SampleIDs = n_distinct(SampleID), .groups = "drop") %>%
  dplyr::filter(n_SampleIDs > 2) %>%
  dplyr::pull(Match_Caseset)
# Step 2: Filter the main data to those problematic casesets,
# keeping relevant columns: SampleID, Match_Caseset, PlateID, Cncr_Caco_Clrt
problematic_data <- data_immuneonc %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::filter(Match_Caseset %in% problematic_casesets) %>%
  dplyr::select(SampleID, Match_Caseset, PlateID, Cncr_Caco_Clrt) %>%
  dplyr::distinct()
# Step 3: Identify Match_Casesets that have multiple PlateIDs (potential problem)
casesets_multiple_plates <- problematic_data %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::summarise(n_PlateIDs = n_distinct(PlateID), .groups = "drop") %>%
  dplyr::filter(n_PlateIDs > 1) %>%
  dplyr::pull(Match_Caseset)
# Step 4: Filter data to only these problematic Match_Casesets with multiple PlateIDs
problematic_multi_plate <- problematic_data %>%
  dplyr::filter(Match_Caseset %in% casesets_multiple_plates)
# Step 5: Find SampleIDs that are the *only* sample on their PlateID within that Match_Caseset (candidates for removal)
removal_candidates <- problematic_multi_plate %>%
  dplyr::group_by(Match_Caseset, PlateID) %>%
  dplyr::filter(n() == 1) %>%
  dplyr::ungroup()
# Step 6: Define a function to check if removing a sample keeps both cases and controls in the Match_Caseset
check_keep_case_control <- function(caseset, sample_to_remove, data) {
  remaining <- data %>%
    dplyr::filter(Match_Caseset == caseset, SampleID != sample_to_remove)
  
  has_case <- any(remaining$Cncr_Caco_Clrt == 1)
  has_control <- any(remaining$Cncr_Caco_Clrt == 0)
  
  return(has_case & has_control)
}
# Step 7: Apply the function to each candidate to test if they can be safely removed
removal_candidates <- removal_candidates %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    can_remove = check_keep_case_control(Match_Caseset, SampleID, problematic_data)
  ) %>%
  dplyr::ungroup()
# Step 8: Filter to samples that can actually be removed safely
removable_samples <- removal_candidates %>%
  dplyr::filter(can_remove)
# Step 9: Output removable samples, or remove them from original data if desired
VAR_sample_exclude <- removable_samples %>%
  dplyr::select(SampleID, Match_Caseset, PlateID)
# Step 10: exclude
data_explore <- data_explore %>%
  dplyr::anti_join(VAR_sample_exclude, by = c("SampleID", "Match_Caseset", "PlateID"))

## exclusion 2: exclude samples missing age or sex ====
data_immuneonc %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::group_by(missing_type = dplyr::case_when(
    is.na(Age_Recr) & is.na(Sex) ~ "Both missing",
    is.na(Age_Recr) ~ "Missing Age_Recr",
    is.na(Sex) ~ "Missing Sex",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!is.na(missing_type)) %>%
  dplyr::distinct(SampleID, .keep_all = TRUE) %>%
  dplyr::count(missing_type)
data_immuneonc <- data_immuneonc %>%
  dplyr::filter(!(Sample_Type == "SAMPLE" & (is.na(Age_Recr) | is.na(Sex))))

data_explore %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::group_by(missing_type = dplyr::case_when(
    is.na(Age_Recr) & is.na(Sex) ~ "Both missing",
    is.na(Age_Recr) ~ "Missing Age_Recr",
    is.na(Sex) ~ "Missing Sex",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!is.na(missing_type)) %>%
  dplyr::distinct(SampleID, .keep_all = TRUE) %>%
  dplyr::count(missing_type)
data_explore <- data_explore %>%
  dplyr::filter(!(Sample_Type == "SAMPLE" & (is.na(Age_Recr) | is.na(Sex))))

## exclusion 3: exclude samples not passing QC ====
data_immuneonc %>% dplyr::select(SampleID, QC_Warning) %>% unique() %>% dplyr::count(QC_Warning)
data_immuneonc <- data_immuneonc %>%
  dplyr::filter(!(Sample_Type == "SAMPLE" & QC_Warning != "Pass"))

data_explore %>% dplyr::select(SampleID, Assay, QC_Warning, Assay_Warning) %>% unique() %>% dplyr::count(QC_Warning, Assay_Warning)
data_explore <- data_explore %>%
  dplyr::filter(!(Sample_Type == "SAMPLE" & (QC_Warning != "PASS" | Assay_Warning != "PASS")))

## exclusion 4: exclude inelligible tumours ====
data_immuneonc %>% dplyr::select(SampleID, Typ_Tumo) %>% unique() %>% dplyr::count(Typ_Tumo)
data_immuneonc <- data_immuneonc %>%
  dplyr::filter(!(Sample_Type == "SAMPLE" & Typ_Tumo %in% c("Exc-Morphology excluded", "Exc-Morphology ineligible")))

data_explore %>% dplyr::select(SampleID, Typ_Tumo) %>% unique() %>% dplyr::count(Typ_Tumo)
data_explore <- data_explore %>%
  dplyr::filter(!(Sample_Type == "SAMPLE" & Typ_Tumo %in% c("Exc-Morphology excluded", "Exc-Morphology ineligible")))

## exclusion 5: for duplicate SampleID with case and control entry, remove the control row ====
sampleids_case_and_control <- data_immuneonc %>%
  filter(Cncr_Caco_Clrt %in% c(0, 1)) %>%
  distinct(SampleID, Cncr_Caco_Clrt) %>%
  count(SampleID) %>%
  filter(n > 1) %>%
  pull(SampleID)
data_immuneonc <- data_immuneonc %>%
  filter(!(SampleID %in% sampleids_case_and_control & Cncr_Caco_Clrt == 0))

sampleids_case_and_control <- data_explore %>%
  filter(Cncr_Caco_Clrt %in% c(0, 1)) %>%
  distinct(SampleID, Cncr_Caco_Clrt) %>%
  count(SampleID) %>%
  filter(n > 1) %>%
  pull(SampleID)
data_explore <- data_explore %>%
  filter(!(SampleID %in% sampleids_case_and_control & Cncr_Caco_Clrt == 0))

## exclusion 6: for duplicate SampleIDs: (1) remove Match_Caseset with only 1 SampleID, (2) remove the Match_Caseset where the duplicate SampleID is a control, (3) remove the first Match_Caseset containing the duplicate SampleID ====
### Step 1: Identify SampleIDs that appear in more than one Match_Caseset
sampleid_multiple_casesets <- dplyr::distinct(data_immuneonc, SampleID, Match_Caseset) %>%
  dplyr::count(SampleID) %>%
  dplyr::filter(n > 1)
### Step 2: Calculate the number of unique SampleIDs per Match_Caseset
caseset_sizes <- dplyr::distinct(data_immuneonc, SampleID, Match_Caseset) %>%
  dplyr::count(Match_Caseset, name = "n_sampleid")
### Step 3: Filter data to only duplicate SampleIDs and keep distinct rows of SampleID, Match_Caseset, and case/control status
dup_info <- data_immuneonc %>%
  dplyr::filter(SampleID %in% sampleid_multiple_casesets$SampleID) %>%
  dplyr::select(SampleID, Match_Caseset, Cncr_Caco_Clrt) %>%
  dplyr::distinct() %>%
  # Join the size of each Match_Caseset (number of unique SampleIDs in it)
  dplyr::left_join(caseset_sizes, by = "Match_Caseset") %>%
  # Focus only on rows where the SampleID is a control (Cncr_Caco_Clrt == 0)
  dplyr::filter(Cncr_Caco_Clrt == 0)
### Step 4: Identify Match_Casesets with only 1 SampleID and remove those Match_Casesets from the main dataset
casesets_to_remove <- dup_info %>%
  dplyr::filter(n_sampleid == 1) %>%
  dplyr::pull(Match_Caseset)
data_immuneonc <- data_immuneonc %>%
  dplyr::filter(!Match_Caseset %in% casesets_to_remove)
### Step 5: Recalculate duplicates after removing small casesets
sampleid_multiple_casesets_after <- data_immuneonc %>%
  dplyr::distinct(SampleID, Match_Caseset) %>%
  dplyr::count(SampleID) %>%
  dplyr::filter(n > 1)
### Step 6: If duplicates remain, exclude :
if (nrow(sampleid_multiple_casesets_after) == 0) {
  # No duplicates remain, finished cleaning
  data_immuneonc <- data_immuneonc
} else {
  # Get info about the remaining duplicates (only controls)
  dup_info_after <- data_immuneonc %>%
    dplyr::filter(SampleID %in% sampleid_multiple_casesets_after$SampleID) %>%
    dplyr::select(SampleID, Match_Caseset, Cncr_Caco_Clrt) %>%
    dplyr::distinct() %>%
    dplyr::left_join(caseset_sizes, by = "Match_Caseset") %>%
    dplyr::filter(Cncr_Caco_Clrt == 0)
  # Identify casesets to remove:
  casesets_to_remove_second <- dup_info_after %>%
    dplyr::group_by(SampleID) %>%
    dplyr::filter(all(n_sampleid == 2)) %>%
    dplyr::slice_min(Match_Caseset, n = 1) %>%
    dplyr::pull(Match_Caseset)
  # Remove those Match_Casesets from the cleaned data
  data_immuneonc <- data_immuneonc %>%
    dplyr::filter(!Match_Caseset %in% casesets_to_remove_second)
}

# Step 1: Identify SampleIDs that appear in more than one Match_Caseset
sampleid_multiple_casesets <- dplyr::distinct(data_explore, SampleID, Match_Caseset) %>%
  dplyr::count(SampleID) %>%
  dplyr::filter(n > 1)
# Step 2: Calculate the number of unique SampleIDs per Match_Caseset
caseset_sizes <- dplyr::distinct(data_explore, SampleID, Match_Caseset) %>%
  dplyr::count(Match_Caseset, name = "n_sampleid")
# Step 3: Filter data to only duplicate SampleIDs and keep distinct rows of SampleID, Match_Caseset, and case/control status
dup_info <- data_explore %>%
  dplyr::filter(SampleID %in% sampleid_multiple_casesets$SampleID) %>%
  dplyr::select(SampleID, Match_Caseset, Cncr_Caco_Clrt) %>%
  dplyr::distinct() %>%
  # Join the size of each Match_Caseset (number of unique SampleIDs in it)
  dplyr::left_join(caseset_sizes, by = "Match_Caseset") %>%
  # Focus only on rows where the SampleID is a control (Cncr_Caco_Clrt == 0)
  dplyr::filter(Cncr_Caco_Clrt == 0)
# Step 4: Identify Match_Casesets with only 1 SampleID and remove those Match_Casesets from the main dataset
casesets_to_remove <- dup_info %>%
  dplyr::filter(n_sampleid == 1) %>%
  dplyr::pull(Match_Caseset)
data_explore <- data_explore %>%
  dplyr::filter(!Match_Caseset %in% casesets_to_remove)
# Step 5: Recalculate duplicates after removing small casesets
sampleid_multiple_casesets_after <- data_explore %>%
  dplyr::distinct(SampleID, Match_Caseset) %>%
  dplyr::count(SampleID) %>%
  dplyr::filter(n > 1)
# Step 6: If duplicates remain, apply second rule:
if (nrow(sampleid_multiple_casesets_after) == 0) {
  # No duplicates remain, finished cleaning
  data_explore <- data_explore
} else {
  # Get info about the remaining duplicates (only controls)
  dup_info_after <- data_explore %>%
    dplyr::filter(SampleID %in% sampleid_multiple_casesets_after$SampleID) %>%
    dplyr::select(SampleID, Match_Caseset, Cncr_Caco_Clrt) %>%
    dplyr::distinct() %>%
    dplyr::left_join(caseset_sizes, by = "Match_Caseset") %>%
    dplyr::filter(Cncr_Caco_Clrt == 0)
  
  # Identify casesets to remove:
  casesets_to_remove_second <- dup_info_after %>%
    dplyr::group_by(SampleID) %>%
    dplyr::filter(all(n_sampleid == 2)) %>%
    dplyr::slice_min(Match_Caseset, n = 1) %>%
    dplyr::pull(Match_Caseset)
  
  # Remove those Match_Casesets from the cleaned data
  data_explore <- data_explore %>%
    dplyr::filter(!Match_Caseset %in% casesets_to_remove_second)
}

## exclusion 7: exclude all Match_Casesets which: (1) do not have at least 2 distinct SampleIDs, (2) that does not have at leas 1 control, and (3) that does not have exactly 1 case ====
check_casesets <- data_immuneonc %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::select(SampleID, Match_Caseset, Cncr_Caco_Clrt) %>%
  unique() %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::summarise(
    n_sampleids = n_distinct(SampleID),
    n_controls = sum(Cncr_Caco_Clrt == 0, na.rm = TRUE),
    n_cases = sum(Cncr_Caco_Clrt == 1, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    passes_criteria = (n_sampleids >= 2) & (n_controls >= 1) & (n_cases == 1)
  )
data_immuneonc <- data_immuneonc %>%
  dplyr::filter(!Match_Caseset %in% (check_casesets %>% dplyr::filter(passes_criteria == FALSE) %>% dplyr::pull(Match_Caseset)))

check_casesets <- data_explore %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::select(SampleID, Match_Caseset, Cncr_Caco_Clrt) %>%
  unique() %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::summarise(
    n_sampleids = n_distinct(SampleID),
    n_controls = sum(Cncr_Caco_Clrt == 0, na.rm = TRUE),
    n_cases = sum(Cncr_Caco_Clrt == 1, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    passes_criteria = (n_sampleids >= 2) & (n_controls >= 1) & (n_cases == 1)
  )
data_explore <- data_explore %>%
  dplyr::filter(!Match_Caseset %in% (check_casesets %>% dplyr::filter(passes_criteria == FALSE) %>% dplyr::pull(Match_Caseset)))

# data checks ====
if ((dup_count <- data_immuneonc %>% select(SampleID, Match_Caseset) %>% unique() %>% count(SampleID) %>% filter(n > 1) %>% nrow()) > 0) {
  message("data_immuneonc: duplicate SampleIDs: ", dup_count)
} else {
  message("data_immuneonc: No duplicate SampleIDs")
}
# Check for Match_Casesets with only 1 unique SampleID
casesets_with_1_sample <- data_immuneonc %>%
  distinct(Match_Caseset, SampleID) %>%
  count(Match_Caseset, name = "n_sampleid") %>%
  filter(n_sampleid == 1)
if (nrow(casesets_with_1_sample) > 0) {
  message("data_immuneonc: Match_Casesets with 1 unique SampleID: ", nrow(casesets_with_1_sample))
} else {
  message("data_immuneonc: Match_Casesets have >1 unique SampleID")
}
# Check that each Match_Caseset has exactly 1 case and 1 control
caseset_case_control_counts <- data_immuneonc %>%
  filter(Sample_Type == "SAMPLE") %>%
  distinct(Match_Caseset, SampleID, Cncr_Caco_Clrt) %>%
  count(Match_Caseset, Cncr_Caco_Clrt) %>%
  tidyr::pivot_wider(names_from = Cncr_Caco_Clrt, values_from = n, values_fill = 0) %>%
  rename(n_control = `0`, n_case = `1`)
casesets_not_1_case_1_control <- caseset_case_control_counts %>%
  filter(n_case != 1 | n_control != 1)
if (nrow(casesets_not_1_case_1_control) > 0) {
  message("data_immuneonc: Match_Casesets do not have exactly 1 case and 1 control: ", nrow(casesets_not_1_case_1_control))
} else {
  message("data_immuneonc: Match_Casesets have exactly 1 case and 1 control")
}




if ((dup_count <- data_explore %>% select(SampleID, Match_Caseset) %>% unique() %>% count(SampleID) %>% filter(n > 1) %>% nrow()) > 0) {
  message("data_explore: duplicate SampleIDs: ", dup_count)
} else {
  message("data_explore: No duplicate SampleIDs")
}
# Check for Match_Casesets with only 1 unique SampleID
casesets_with_1_sample <- data_explore %>%
  distinct(Match_Caseset, SampleID) %>%
  count(Match_Caseset, name = "n_sampleid") %>%
  filter(n_sampleid == 1)
if (nrow(casesets_with_1_sample) > 0) {
  message("data_explore: Match_Casesets with 1 unique SampleID: ", nrow(casesets_with_1_sample))
} else {
  message("data_explore: Match_Casesets have >1 unique SampleID")
}
# Check that each Match_Caseset has exactly 1 case and 1 control
caseset_case_control_counts <- data_explore %>%
  filter(Sample_Type == "SAMPLE") %>%
  distinct(Match_Caseset, SampleID, Cncr_Caco_Clrt) %>%
  count(Match_Caseset, Cncr_Caco_Clrt) %>%
  tidyr::pivot_wider(names_from = Cncr_Caco_Clrt, values_from = n, values_fill = 0) %>%
  rename(n_control = `0`, n_case = `1`)
casesets_not_1_case_1_control <- caseset_case_control_counts %>%
  filter(n_case != 1 | n_control != 1)
if (nrow(casesets_not_1_case_1_control) > 0) {
  message("data_explore: Match_Casesets do not have exactly 1 case and 1 control: ", nrow(casesets_not_1_case_1_control))
} else {
  message("data_explore: Match_Casesets have exactly 1 case and 1 control")
}

## write ====
write.table(data_immuneonc, "data/processed/EPIC_olink-immuneonc.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(data_explore, "data/processed/EPIC_olink-explore.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

