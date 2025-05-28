#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------------------------------------------------
# R script : Preprocessing step to obtain Robject to continue the analysis
# Auteur  : Léa ROGUE
# Date    : 23-01-2025
# Description : This script use oligo to process microarrays data. Firstly, the script read all the files and create the 
# object with metadata from the .sdrf file to have all the informations. Finally, object of expression from raw data are
# saved to purchase quality analysis with graphics for exemple (pca or boxplot). The expression data are normalized with 
# the RMA algorithm (from oligo) to delete the background noise and normalize data, here an object is also saved.
# ----------------------------------------------------------------------------------------------------------------------

# If an error occur for the rma fonction, these commands line solve it
#BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
#BiocManager::install("oligo", configure.args="--disable-threading", force = TRUE)

# Load packages
library(oligo)

# Directories paths
data_dir <- "../data/E-MTAB-7264_full"
output_dir <- "../results"
dir.create(output_dir)

# SRDF to store the metadata
SDRF <- read.delim(paste0(data_dir, "/E-MTAB-7264.sdrf.txt"))

# Give raw names the sample name 
rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)

# Read the raw data and associate metadata
raw_data <- read.celfiles(filenames = file.path(data_dir, SDRF$Array.Data.File), phenoData = SDRF)
message("Raw data object created")

# Save raw data object
message("Object raw data saving...")
save(raw_data, file = paste0(output_dir, "/object_chondrosarcoma_raw.RData"))
message("Object saved")

# Normalization and delete background noise
norm_data <- rma(raw_data, target = "core")

# Save normalized data
message('Object saving...')
save(norm_data, file = paste0(output_dir, "/object_chondrosarcoma_norm.RData"))
message('Object saved')
