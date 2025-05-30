#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------------------------------------------------
# R script : Data processing of gene expression in normal tissues from GTEx database of CTA 
# Author  : Léa ROGUE
# Date    : 24-03-2025
# Description : This scripts is to process GTEx data in order to obtain mean of the expression for each CTA in each normal
# tissue. This allow us to see the expression of CTAs of interest in normal tissues to select the tumor specific CTAs.
# ------------------------------------------------------------------------------------------------------------------------

# Download necessary files with bash commands
#wget https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt
#wget https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz

# Load libraries
library(dplyr)
library(data.table)

# Read samples attributes
df_samples <- read.delim("../data/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt", as.is=TRUE, header=TRUE)

# Select RNAseq samples
df_rnaseq <- df_samples[grep("GTEX", df_samples$SAMPID) & df_samples$SMAFRZE == "RNASEQ", c("SAMPID", "SMTS")]

# Read TPM
df_tpm <- fread("../data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz")

# Select CTA
cta <- (as.vector(read.table("../data/CTA_list_clean.txt", header = FALSE)))$V1
housekp_genes <- c("EGFR", "ERBB2", "FOLH1", "SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5")

# Select data in TPM matrix
df_expr_cta <- df_tpm[df_tpm$Description %in% cta | df_tpm$Description %in% housekp_genes,]
cta_names <- df_expr_cta$Description
df_expr_cta <- as.data.frame(df_expr_cta)
df_expr_cta <- df_expr_cta[, colnames(df_expr_cta) %in% df_rnaseq$SAMPID]

# Convert and merge the data with tissue type
df_expr_long <- as.data.frame(t(df_expr_cta))
colnames(df_expr_long) <- cta_names
df_expr_long$SAMPID <- rownames(df_expr_long)
df_expr_cta_tissues <- merge(df_expr_long, df_rnaseq, by = "SAMPID")

# Compute mean for each tissue
df_expr_cta_tissues <- df_expr_cta_tissues[, -1]
df_tissue_cta_mean <- df_expr_cta_tissues %>%
  group_by(SMTS) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) 

# Save
write.table(df_tissue_cta_mean, "../results/matrix_expr_cta_tissues_gtex.tsv",  sep = "\t", row.names = FALSE, quote = FALSE)
