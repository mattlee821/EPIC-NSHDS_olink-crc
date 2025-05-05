# EPIC-NSHDS_olink-CRC

This project investigates the association between proteins measured using the Olink platform and colorectal cancer in the European Prospective Investigation into Cancer and Nutrition study (EPIC) and the Northern Sweden Health and Disease Study (NSHDS)

---

# To do
1. heterogeneity analysis
2. EPIC explore analysis
2. NSHDS analysis
3. paper: re-format and change structure for NSHDS inclusion
  - figures
  - tables

# Done
1. EPIC immuneonc analysis
2. MR/colocalisation

# Project  Structure  
This project is organised into the following directories:  

- `data/`
    - contains raw and processed data
    - `src/001_data-processing/001_data-process.R` takes input from `raw/` and outputs into `processed/`
- `src/`
    - contains all code
    - all code output (except as stated above for `data/`) is stored in `analysis/` in directories of the same name
    - use `000_source.R` to install packages
    - `000_funcistion.R` contains all custom functions used
- `analysis/`
    - contains the output from `src/` in directories of the same name

You can make the directory structure by copying the file `directory-structure.txt` and running the below code:
```
# Recreate the directory structure in target_dir
directories <- data.table::fread("directory-structure.txt")
target_dir <- "EPIC-NSHDS_olink-crc/"
# Ensure target directory exists
if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}
# Recreate the directory structure under target_dir
for (i in directories$dir) {
  new_path <- file.path(target_dir, i)
  dir.create(new_path, recursive = TRUE, showWarnings = FALSE)
}
```

- TBC
    - `code_review/` – Helps a "code buddy" reproduce manuscripts, presentations, and reports.
    - `manuscripts/`, `presentations/`, and `reports/` – Store outputs for different sub-projects.

#### `001_data-processing/`
- provides scripts for initial processing of raw data
- `001_data-process.R` loads raw data, performs some general exclusions, and formats and creates dataframes to be used in subsequent scripts
    - the output of this should be 3 data frames: (1) feature data, which contains column 1 as IDs and all remaining columns as features; feature meta-data, which contains columns with the names of all features and any associated feature meta-data (e.g., LOD); (3) sample meta-data, which contains column 1 as IDs and all remaining columns as meta-data (e.g., age, sex, cancer status)
- `002_data-check*.R` loads the output of `001_data-process.R` and calculates coefficient of variation, inter-class correlations, and plots distributions and scatter plots for every feature
- `003_make-processed-data.R` performs omics data processing (`OmicsProcessing::process_data()`)
    - this script needs to be run twice, once using the data loaded on line 14 and once with the df created from line 77
    - you then need to change the arguments for `imputation` and `transformation` to include all combinations of `imputation = FALSE/LOD` and `transformation = FALSE/InvRank` for both dataframes above
    - the output of this script will be 1 `.rda` file for each combination of the values above (4 `.rda` files for each data frame, 8 files in total)

#### `002_EPIC-analysis/`
- `001_analysis/001_master.R` loads the output from `003_make-processed-data.R` and the associated meta-data and then runs the regression analyses
- XXX
- XXX

