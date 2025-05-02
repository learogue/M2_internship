#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------
# Python script : Analyze the number of peptides per gene
# Author  : Léa ROGUE
# Date    : 02-05-2025
# Description : This script analyzes the number of peptides per gene from the results of the strong binding peptides for
# HLA-A0201.
# ----------------------------------------------------------------------------------------------------------------------

from collections import defaultdict 
import pandas as pd

# Read gene names and their corresponding entry names
d_pep_gene = defaultdict(list)
d_genes = defaultdict()
with open('../data/entry_name_gene_name.tsv') as f:
    for lig in f:
        lig = lig.strip().split()
        d_genes[lig[0]] = lig[1]

# Read peptides IDs and their corresponding gene
with open('../results/selected_results_blastp_mismatch.tsv') as f:
    for lig in f:
        lig = lig.strip().split()
        d_pep_gene[lig[0]].append(d_genes[lig[4]])

# Count the number of peptides per gene
d_nb_pep_gene = defaultdict(int)
with open('../results/strong_binding_peptides_hla_a0201.out') as f:
    for lig in f:
        d_nb_pep_gene[str(d_pep_gene[lig.strip().split()[10]]).replace('\'', '').replace('[', '').replace(']', '')] += 1

# Sorted dictionary
sorted_dict = dict(sorted(d_nb_pep_gene.items(), key=lambda item: item[1], reverse=True))

# Create dataframe
df = pd.DataFrame(list(sorted_dict.items()), columns=['Genes', 'Peptides number'])

# Save
df.to_csv('../results/nb_pep_per_gene_sb_hla_a0201.tsv', index=False, sep = '\t')
