library(data.table)

# Define source directory
source_dir <- "/data/IARC_Biostat/work/EPIC-NSHDS_olink-crc/"

# Get all subdirectories (excluding the root itself)
dir_paths <- list.dirs(path = source_dir, full.names = TRUE, recursive = TRUE)
dir_paths <- dir_paths[dir_paths != source_dir]

# Exclude hidden directories and '/code-review'
dir_paths <- dir_paths[!grepl("/\\.", dir_paths) & !grepl("/code_review", dir_paths)]

# Convert to relative paths
relative_paths <- gsub(pattern = source_dir, replacement = "", x = dir_paths, fixed = TRUE)

# Remove any empty strings
relative_paths <- relative_paths[relative_paths != ""]

# Save to file
dir_df <- data.frame(dir = relative_paths, stringsAsFactors = FALSE)
write.table(dir_df, file = "directory-structure.txt", sep = "\t", col.names = T, row.names = F)

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

