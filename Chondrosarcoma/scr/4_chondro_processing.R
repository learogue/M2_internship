#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------------------------------------------------
# R script : Data processing of gene expression
# Auteur  : Léa ROGUE
# Date    : 12-02-2025
# Description : This script processes and analyzes normalized gene expression data. It integrates immune cell signature 
# genes and CTA (Cancer-Testis Antigen) genes and calculates Z-scores.
# Steps:
#   1. Load the normalized expression data.
#   2. Rename patient ID columns and associate probe IDs with gene symbols.
#   3. Merge expression data with CTA and immune cell signature gene lists.
#   4. Clean and organize the dataset by consolidating duplicate probes and calculating mean expression values.
#   5. Compute Z-scores to standardize gene expression levels.
# ------------------------------------------------------------------------------------------------------------------------

# Load packages
library(oligo)
library(hugene20sttranscriptcluster.db)
library(dplyr)

# Function to calculate Z scores
calculate_z_scores <- function(df_input, col) {
  # Exclude col PROBEID SYMBOL and CTA
  data_values <- df_input[, -c(col)]

  # Calculate Z-scores
  z_scores_row <- t(scale(t(data_values)))

  # Add columns
  df_z_scores <- cbind(df_input[, c(col)], z_scores_row)

  # Return df
  return(df_z_scores)
}

# Load the normalized data
load("../results/object_chondrosarcoma_norm.RData")

# Save the metadata
metadata_col <- c("Characteristics.exp.subtypes.",
	"Characteristics.mom.subtypes.",
	"Characteristics.idh1.aamut.",
	"Characteristics.idh1.freq.",
	"Characteristics.idh2.aamut.",
	"Characteristics.idh2.freq.",
	"Characteristics.col2a1.",
	"Characteristics.tp53.",
	"Characteristics.cdkn2a.copynumber.",
	"Characteristics.os.delay.",
	"Unit.time.unit.",
	"Characteristics.event.death.",
	"Characteristics.tumor.grading.",
	"Factor.Value.tumor.grading.")
df <- pData(norm_data)[,metadata_col]
write.table(df, file = "../results/metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Take only the coding genes
mapped_probes <- mappedkeys(hugene20sttranscriptclusterGENENAME)

# Associate genes with probes
anno <- AnnotationDbi::select(hugene20sttranscriptcluster.db, keys = mapped_probes, keytype = "PROBEID", columns = c("SYMBOL", "GENENAME"))

# Take rows that have a gene
anno <- subset(anno, !is.na(SYMBOL))

# Change names of column with patient IDs
colnames(norm_data) <- norm_data$Source.Name
rownames(pData(norm_data)) <- norm_data$Source.Name

# Expression data and rename the rows with probe id
expr_norm_data_df <- as.data.frame(exprs(norm_data))
expr_norm_data_df$PROBEID <- rownames(expr_norm_data_df)

# Merge expression data with SYMBOL
merged_df <- merge(anno, expr_norm_data_df, by = "PROBEID", all.x = TRUE)
merged_df <- merged_df[, c(-1,-3)]

# For genes that have multiple probe, calculate the mean
merged_df_avg <- merged_df %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# Read CTA file
merged_df_avg$CTA <-  NA
df_CTA <- read.table("../data/CTA_list_clean.txt", header = FALSE)

# Rename col
colnames(df_CTA) <- c("SYMBOL")

# Df with CTA and whole genes by adding CTA in col CTA for CTA genes
df_CTA_whole <- merged_df_avg %>%
  left_join(df_CTA, by = "SYMBOL") %>%
  mutate(CTA = ifelse(SYMBOL %in% df_CTA$SYMBOL, "CTA", NA))

# Reorganize columns
df_CTA_whole <- df_CTA_whole %>%
  select(SYMBOL, CTA, everything())

# Read immune cells genes
df_immune_sign <- read.table("../data/immune_cells_genes.tsv", header = FALSE, sep = "\t")
colnames(df_immune_sign) <- c("Signature", "Gene")

# Merge df_CTA_whole and df_immune_sign to associate signatures
df_CTA_immune_sign_whole <- merge(df_CTA_whole, df_immune_sign, by.x = "SYMBOL", by.y = "Gene", all.x = TRUE)

# Reorganize the columns if needed
df_CTA_immune_sign_whole <- df_CTA_immune_sign_whole %>%
  select(SYMBOL, CTA, Signature, everything())

# Combine multiple signatures into one for each gene by concatenating with a comma
df_CTA_immune_sign_whole_clean <- df_CTA_immune_sign_whole %>%
  group_by(SYMBOL) %>%
  summarise(
    Signature = paste(unique(Signature), collapse = ", "),                                              
    across(everything(), ~first(.)),  
    .groups = "drop")
write.table(df_CTA_immune_sign_whole_clean, file = "../results/whole_gene_int_CTA_sign_imm_clean.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Z-scores on rows on whole genes to comparate expression between genes with a robust manner
df_CTA_immune_whole_clean_z_scores <- calculate_z_scores(df_CTA_immune_sign_whole_clean, c(1,2,3))
write.table(df_CTA_immune_whole_clean_z_scores, file = "../results/whole_gene_CTA_sign_imm_clean_avg_z_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Z-scores on rows on whole genes to comparate expression between genes with a robust manner
df_CTA_immune_whole_z_scores <- calculate_z_scores(df_CTA_immune_sign_whole, c(1,2,3))
write.table(df_CTA_immune_whole_z_scores, file = "../results/whole_gene_CTA_sign_imm_avg_z_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Average the expression between same immune cells types
# Take rows with immune cells signature from normalized data
df_avg_immune_sign <- df_CTA_immune_sign_whole %>%
  filter(Signature != "NA")

# Group by signature and calculate mean of expression values
df_avg_immune_sign_final <- df_avg_immune_sign %>%
  select(-c(SYMBOL, CTA)) %>% 
  group_by(Signature) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) 

# Z-scores
df_avg_immune_sign_z_scores <- calculate_z_scores(df_avg_immune_sign_final, 1)
write.table(df_avg_immune_sign_z_scores, file = "../results/imm_sign_avg_z_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Read MHC genes file
df_MHC <- read.table("../data/MHC_genes.txt", header = TRUE, sep = "\t")

# Merge the df
df_expr_MHC <- merge(df_MHC, df_CTA_immune_sign_whole_clean, by.x = "SYMBOL")

# Average expression for same MHC type
df_MHC <- df_expr_MHC %>%
  select(-c(SYMBOL, Signature, CTA)) %>%
  group_by(Type) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
write.table(df_expr_MHC, file = "../results/expr_MHC.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Z scores on MHC genes
df_expr_MHC_z_scores <- calculate_z_scores(df_expr_MHC, c(1,2,3,4))
write.table(df_expr_MHC_z_scores, file = "../results/expr_MHC_z_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Immune cells and MHC genes
# Read immune cells genes
df_immune_sign <- read.table("../data/immune_cells_mhc_genes.tsv", header = FALSE, sep = "\t")
colnames(df_immune_sign) <- c("Signature", "Gene")

# Merge df_CTA_whole and df_immune_sign to associate signatures
df_CTA_immune_sign_whole <- merge(df_CTA_whole, df_immune_sign, by.x = "SYMBOL", by.y = "Gene", all.x = TRUE)

# Reorganize the columns if needed
df_CTA_immune_sign_whole <- df_CTA_immune_sign_whole %>%
  select(SYMBOL, CTA, Signature, everything())

# Average the expression between same immune cells types
# Take rows with immune cells signature from normalized data
df_avg_immune_sign <- df_CTA_immune_sign_whole %>%
  filter(Signature != "NA")

# Group by signature and calculate mean of expression values
df_avg_immune_sign_final <- df_avg_immune_sign %>%
  select(-c(SYMBOL, CTA)) %>% 
  group_by(Signature) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) 

# Z-scores
df_avg_immune_sign_z_scores <- calculate_z_scores(df_avg_immune_sign_final, 1)
write.table(df_avg_immune_sign_z_scores, file = "../results/imm_sign_mhc_genes_avg_z_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
