#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------------------------
# R script : Compute transcripts per million from pseudo-bulk scRNAseq
# Auteur  : Léa ROGUE
# Date    : 26-03-2025
# Description : This python script read anndata object, extract counts and summary it per samples. Next, lengths of genes are search and
# calculated in gtf file. Then TPM are computed and saved
# ----------------------------------------------------------------------------------------------------------------------------------------

# Load libraries
import pandas as pd
import scanpy as sc
import anndata as ad
from collections import defaultdict
import re
from bioinfokit.analys import norm

# Read object created with scRNAseq tool from my M1 internship
adata = ad.read_h5ad("../data/object_merged_2.h5ad")

# Aggregate counts
aggr = sc.get.aggregate(adata, by = 'dataset', func = 'sum')

# Transform and save the df
df = aggr.to_df(layer="sum")
#df.to_csv('matrix_pseudo_bulk.tsv', sep='\t', index=True)

# List of genes
l = list(df.columns)

# Read GTF file
d_length = defaultdict(int)
with open('../data/Homo_sapiens.GRCh38.113.gtf', 'r') as f:
    for lig in f:
        lig = lig.strip()
        lig = lig.split('\t')
        if len(lig) == 9 and lig[2] == 'exon':
            match = re.search(r'gene_name "([^"]+)"', lig[8])
            if match:
                gene_name = match.group(1)  # Extract name
                if gene_name in l:
                    d_length[gene_name] += abs(int(lig[4]) - int(lig[3]) + 1)

# Convert
lengths = pd.Series(d_length)

# Transposer df pour que les gènes soient en index
df_transposed = df.T

# Filter
df_transposed = df_transposed[df_transposed.index.isin(d_length.keys())]

# Add length column and save
df_transposed['length'] = df_transposed.index.map(lengths)
df_transposed.to_csv("../results/matrix_pseudo_bulk_length.tsv", sep='\t', index=True)

# TPM computing and save
nm = norm()
tpm_norm = nm.tpm(df=df_transposed, gl='length')
df_tpm = nm.tpm_norm
df_tpm.to_csv("../results/matrix_pseudo_bulk_tpm_normalized.tsv", sep='\t', index=True)
