rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)

# data_samples ====
data_samples <- haven::read_sas(data_file = "data/raw/EPIC_phenotype/clrt_caco.sas7bdat")

# EPIC olink-immuneonc ====
## proteomics sample meta-data ====

## VARS 
VAR_control <- "CONTROL"
VAR_id <- "Idepic_Bio"

## data 
data <- xlsx::read.xlsx(file = "data/raw/EPIC_olink-immuneonc/20211440_Murphy_NPX.xlsx", sheetIndex = 1)

## select just meta data
data_meta <- data %>%
  dplyr::slice(-c(1:6, 2052:nrow(data))) %>%
  dplyr::select(1,94,95) %>%
  tibble::as_tibble()
names(data_meta) <- c("Idepic_Bio", "batch_plate", "qc")

## select controls    
data_meta_controls <- data_meta %>%
  dplyr::filter(if_any(everything(), ~ grepl(VAR_control, .)))

## join phenofile
data_meta <- data_meta %>%
  dplyr::filter(!if_any(everything(), ~ grepl(VAR_control, .))) %>%
  dplyr::left_join(data_samples, by = VAR_id, relationship = "many-to-many") 

## exclusions
### step 1: exclude samples not passing QC
data_meta %>% dplyr::count(qc)
data_meta <- data_meta %>%
  dplyr::filter(qc == "Pass") 

### step 2: exclude inelligible tumours
data_meta %>% dplyr::count(Typ_Tumo)
data_meta <- data_meta %>%
  dplyr::filter(!(Typ_Tumo %in% c("Exc-Morphology excluded", "Exc-Morphology ineligible")))

### step 3: exclude individuals without a matched caseset
data_meta %>%
  dplyr::count(Match_Caseset) %>%
  dplyr::filter(n == 1) %>%
  dplyr::summarise(IDs = n())

data_meta <- data_meta %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::filter(n() >= 2) %>%
  dplyr::ungroup()

### Step 4: exclude matched casesets with more than two people
data_meta %>%
  dplyr::count(Match_Caseset) %>%
  dplyr::filter(n > 2) %>%
  dplyr::summarise(IDs = n())

data_meta <- data_meta %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::filter(n() <= 2) %>%
  dplyr::ungroup()

### step 5: exclude IDs present in more than one Match_Caseset 
data_meta %>%
  dplyr::distinct(Idepic, Match_Caseset) %>%
  dplyr::count(Idepic) %>%       
  dplyr::filter(n > 1) %>%       
  dplyr::summarise(IDs = n())        

data_meta <- data_meta %>%
  dplyr::filter(!Idepic %in% (data_meta %>%
                                dplyr::group_by(Idepic) %>%
                                dplyr::summarise(unique_casesets = dplyr::n_distinct(Match_Caseset)) %>%
                                dplyr::filter(unique_casesets > 1) %>%
                                dplyr::pull(Idepic)))

### step 6: exclude individuals without a matched caseset
data_meta %>%
  dplyr::count(Match_Caseset) %>%
  dplyr::filter(n == 1) %>%
  dplyr::summarise(IDs = n())

data_meta <- data_meta %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::filter(n() >= 2) %>%
  dplyr::ungroup()

## write
write.table(data_meta, "data/processed/EPIC_olink-immuneonc_sample-metadata.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(data_meta_controls, "data/processed/EPIC_olink-immuneonc_sample-metadata-controls.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## proteomics feature meta-data ====

## VARS 
VAR_control <- "CONTROL"
VAR_id <- "Idepic_Bio"

## data 
data <- readxl::read_xlsx("data/raw/EPIC_olink-immuneonc/20211440_Murphy_NPX.xlsx")

## select just meta data
data_meta <- data %>%
  dplyr::slice(c(3:5, 2054:nrow(data))) %>% 
  t() %>%
  tibble::as_tibble() %>%
  dplyr::filter(!if_any(everything(), ~ grepl("Plate ID|QC Warning", .))) 
colnames(data_meta) <- data_meta[1, ] # Set the first row as column names
data_meta <- data_meta[-1, ] %>% # Remove the first row and reset row names
  tibble::rownames_to_column(var = "rownames") %>%
  dplyr::select(-rownames) %>%
  dplyr::rename(UNIPROT = `Uniprot ID`,
         missing_pct = `Missing Data freq.`) %>%
  dplyr::mutate(missing_pct = as.numeric(missing_pct))

## write 
write.table(data_meta, "data/processed/EPIC_olink-immuneonc_feature-metadata.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## data proteomics samples-features ====

## VARS 
VAR_control <- "CONTROL"
VAR_id <- "Idepic_Bio"

## data 
data <- readxl::read_xlsx("data/raw/EPIC_olink-immuneonc/20211440_Murphy_NPX.xlsx")
data_meta_samples <- functions::read_file(file_path = "data/processed/EPIC_olink-immuneonc_sample-metadata.txt")

## select data 
data <- data %>%
  slice(-c(1:4, 6, 2053:nrow(data))) %>%
  select(-c(94, 95)) %>%
  as_tibble() 
colnames(data) <- data[1, ] 
data <- data[-1, ] %>%
  rename(!!VAR_id := names(data)[1])
data_controls <- data %>%
  filter(if_any(everything(), ~ grepl(VAR_control, .)))
data <- data %>%
  filter(Idepic_Bio %in% data_meta_samples$Idepic_Bio)

## write
write.table(data, "data/processed/EPIC_olink-immuneonc_feature-data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(data_controls, "data/processed/EPIC_olink-immuneonc_feature-data-controls.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



# EPIC olink-explore ====
## proteomics sample meta-data ====

## VARS 
VAR_id <- "Idepic_Bio"

## data 
data <- read.table(file = "data/raw/EPIC_olink-explore/WL-3904_NPX_2024-12-09.csv", sep = ";", header = T)

## select just meta data
data_meta <- data %>%
  dplyr::select(SampleID, PlateID, Sample_Type) %>%
  tibble::as_tibble()
names(data_meta) <- c("Idepic_Bio", "batch_plate", "Sample_Type")

## select controls    
data_meta_controls <- data_meta %>%
  dplyr::filter(Sample_Type == "CONTROL")

## join phenofile
data_meta <- data_meta %>%
  dplyr::filter(Sample_Type == "SAMPLE") %>%
  dplyr::select(-Sample_Type) %>%
  unique() %>%
  dplyr::left_join(data_samples, by = VAR_id, relationship = "many-to-many") 

## exclusions
### step 1: exclude inelligible tumours
data_meta %>% dplyr::count(Typ_Tumo)
data_meta <- data_meta %>%
  dplyr::filter(!(Typ_Tumo %in% c("Exc-Morphology excluded", "Exc-Morphology ineligible")))

### step 2: exclude individuals without a matched caseset
data_meta %>%
  dplyr::count(Match_Caseset) %>%
  dplyr::filter(n == 1) %>%
  dplyr::summarise(IDs = n())

data_meta <- data_meta %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::filter(n() >= 2) %>%
  dplyr::ungroup()

### Step 3: exclude matched casesets with more than two people
data_meta %>%
  dplyr::count(Match_Caseset) %>%
  dplyr::filter(n > 2) %>%
  dplyr::summarise(IDs = n())

data_meta <- data_meta %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::filter(n() <= 2) %>%
  dplyr::ungroup()

### step 4: exclude IDs present in more than one Match_Caseset 
data_meta %>%
  dplyr::distinct(Idepic, Match_Caseset) %>%
  dplyr::count(Idepic) %>%       
  dplyr::filter(n > 1) %>%       
  dplyr::summarise(IDs = n())        

data_meta <- data_meta %>%
  dplyr::filter(!Idepic %in% (data_meta %>%
                                dplyr::group_by(Idepic) %>%
                                dplyr::summarise(unique_casesets = dplyr::n_distinct(Match_Caseset)) %>%
                                dplyr::filter(unique_casesets > 1) %>%
                                dplyr::pull(Idepic)))

### step 5: exclude individuals without a matched caseset
data_meta %>%
  dplyr::count(Match_Caseset) %>%
  dplyr::filter(n == 1) %>%
  dplyr::summarise(IDs = n())

data_meta <- data_meta %>%
  dplyr::group_by(Match_Caseset) %>%
  dplyr::filter(n() >= 2) %>%
  dplyr::ungroup()

## write
write.table(data_meta, "data/processed/EPIC_olink-explore_sample-metadata.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(data_meta_controls, "data/processed/EPIC_olink-explore_sample-metadata-controls.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## proteomics feature meta-data ====

## data 
data <- read.table(file = "data/raw/EPIC_olink-explore/WL-3904_NPX_2024-12-09.csv", sep = ";", header = T)

## select just meta data - there is an LOD per plate
data_meta <- data %>%
  dplyr::select(OlinkID, UniProt, Assay, MissingFreq, LOD, PlateID, Panel, Assay_Warning) %>%
  unique() %>%
  tibble::as_tibble() %>%
  dplyr::rename(missing_pct = `MissingFreq`,
                batch_plate = PlateID,
                qc = Assay_Warning) %>%
  dplyr::mutate(missing_pct = as.numeric(missing_pct))

## write 
write.table(data_meta, "data/processed/EPIC_olink-explore_feature-metadata.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## data proteomics samples-features ====

## VARS 
VAR_id <- "Idepic_Bio"

## data 
data <- read.table(file = "data/raw/EPIC_olink-explore/WL-3904_NPX_2024-12-09.csv", sep = ";", header = T)
data_meta_samples <- functions::read_file(file_path = "data/processed/EPIC_olink-explore_sample-metadata.txt")

## select data 
data <- data %>%
  dplyr::select(SampleID, OlinkID, NPX) %>%
  tibble::as_tibble() %>%
  tidyr::pivot_wider(names_from = OlinkID, values_from = NPX)

data_controls <- data %>%
  dplyr::filter(if_any(everything(), ~ grepl("Run", .)))
data <- data %>%
  dplyr::filter(!if_any(everything(), ~ grepl("Run", .)))

## write
write.table(data, "data/processed/EPIC_olink-explore_feature-data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(data_controls, "data/processed/EPIC_olink-explore_feature-data-controls.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


