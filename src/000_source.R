# packages: conditional install ====
## Load remotes for GitHub installations
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

## Load BiocManager for Bioconductor installations
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## Define CRAN, GitHub, and Bioconductor packages separately
cran_packages <- c("data.table", "dplyr", "tidyr", "survival", "ggplot2", "haven")
github_packages <- c("NightingaleHealth/ggforestplot", "IARCBiostat/OmicsProcessing", "IARCBiostat/ImputationReport", "privefl/bigutilsr")
bioc_packages <- c("pcaMethods", "impute", "imputeLCMD")

## Install missing CRAN packages
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(new_cran_packages) > 0) {
  install.packages(new_cran_packages)
}

## Install missing GitHub packages
for (pkg in github_packages) {
  pkg_name <- basename(pkg)  # Extract package name from repo path
  if (!pkg_name %in% installed.packages()[, "Package"]) {
    remotes::install_github(pkg)
  }
}

## Install missing Bioconductor packages
new_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if (length(new_bioc_packages) > 0) {
  BiocManager::install(new_bioc_packages)
}

# libraries: load ====
lapply(c(cran_packages, "plinkbinr"), library, character.only = TRUE)

