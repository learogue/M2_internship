#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------------------------------------------------
# R script : Script to use ArrayQualityMetrics to perform quality control of the data
# Auteur  : Léa ROGUE
# Date    : 22-01-2025
# Description : This script perform quality control. It needs a lot of informatic resources because of 102 files. This 
# generates figures and a html page to summarize the quality control.
# ------------------------------------------------------------------------------------------------------------------------

# Load packages
library(oligo)
library(arrayQualityMetrics)

# Loading the objects
load("../results/object_chondrosarcoma_raw.RData")
load("../results/object_chondrosarcoma_norm.RData")
message("Object loaded")

# Array quality control
arrayQualityMetrics(expressionset = raw_data, outdir = "../results/qc/qc_raw", force = TRUE, do.logtransform = TRUE)
arrayQualityMetrics(expressionset = norm_data, outdir = "../results/qc/qc_norm", force = TRUE, do.logtransform = TRUE)
message("Quality report finished")
