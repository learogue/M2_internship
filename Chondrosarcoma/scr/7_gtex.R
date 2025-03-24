

#wget https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt
#wget https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz

library(DESeq2)
library(dplyr)
library(data.table)

# Read samples attributes
df_samples <- read.delim("../data/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt", as.is=TRUE, header=TRUE)

# Select RNAseq samples
df_rnaseq <- df_samples[grep("GTEX", df_samples$SAMPID) & df_samples$SMAFRZE == "RNASEQ" & df_samples$SMTS != "Testis",
  c("SAMPID", "SMTS")]

# Read tpm
df_tpm <- fread("../data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz")

# Select CTA
cta <- (as.vector(read.table("../data/CTA_list_clean.txt", header = FALSE)))$V1

df_expr_cta <- df_tpm[df_tpm$Description %in% cta,]
cta_names <- df_expr_cta$Description

df_expr_cta <- as.data.frame(df_expr_cta)
df_expr_cta <- df_expr_cta[, colnames(df_expr_cta) %in% df_rnaseq$SAMPID]

# Convertir la matrice en un format "long" pour faciliter l'agrégation
df_expr_long <- as.data.frame(t(df_expr_cta))
colnames(df_expr_long) <- cta_names
df_expr_long$SAMPID <- rownames(df_expr_long)
df_expr_cta_tissues <- merge(df_expr_long, df_rnaseq, by = "SAMPID")

# Calculer la moyenne de l'expression pour chaque tissu
df_expr_cta_tissues <- df_expr_cta_tissues[, -1]
df_tissue_cta_mean <- df_expr_cta_tissues %>%
  group_by(SMTS) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) 

#write.table(df_tissue_cta_mean, "../matrix_expr_cta_tissues_gtex.tsv",  sep = "\t", row.names = FALSE, quote = FALSE)

