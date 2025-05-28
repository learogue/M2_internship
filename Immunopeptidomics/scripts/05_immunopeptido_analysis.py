#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------
# Python script : Analyze the affinity per gene and integrate it with the expression level
# Author  : Léa ROGUE
# Date    : 02-05-2025
# Description : This script analyzes the affinity of peptides per gene from the results of the strong binding peptides for
# HLA-A0201.
# ----------------------------------------------------------------------------------------------------------------------

# Import packages
from collections import defaultdict 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import plotly.graph_objects as go

# Read gene names and their corresponding entry names
d_pep_gene = defaultdict(list)
d_genes = defaultdict()
with open('../data/entry_name_gene_name.tsv', 'r') as f:
    for lig in f:
        lig = lig.strip().split()
        d_genes[lig[0]] = lig[1]

# Read peptides IDs and their corresponding gene
with open('../results/selected_results_blastp_0_mismatch.tsv', 'r') as f:
    for lig in f:
        lig = lig.strip().split()
        d_pep_gene[lig[0]].append(d_genes[lig[4]])

# Adding column with genes in the file with genes
with open('../results/netmhc_sb/res_netmhc_HLA-A0201_selected_blastp_sb_annotated.out', 'w') as f_w:
    with open('../results/netmhc_sb/res_netmhc_HLA-A0201_selected_blastp_sb.out', 'r') as f_r:
        for lig in f_r:
            genes = str(d_pep_gene[lig.strip().split()[10]]).replace('\'', '').replace('[', '').replace(']', '')
            lig = lig.strip() + '\t' + genes
            f_w.write(lig + '\n')

# Read the data
df_int = pd.read_csv('../../Chondrosarcoma/results/whole_gene_int_CTA_sign_imm_clean.tsv', sep = '\t')

# Select CTA genes
df_int = df_int[df_int['SYMBOL'].isin(d_genes.values())]

# Mean expression level
df_int['Mean_expression'] = df_int.iloc[:,3:].mean(axis = 1)

# Have affinity for each peptide
# Prepare data
d_pep_affinity = defaultdict(dict)
with open('../results/netmhc_sb/res_netmhc_HLA-A0201_selected_blastp_sb_annotated.out', 'r') as f:
    for lig in f:
        lig = lig.strip().split()
        d_pep_affinity[lig[2]][' '.join(lig[16:])] = float(lig[12])

# Convert to DataFrame with 3 columns: Gene, Peptide, Affinity
l = []
for pep, d_pep in d_pep_affinity.items():
    for gene, affinity in d_pep.items():
        l.append((pep, gene, affinity))

# Create and DataFrame
df_pep_genes_aff = pd.DataFrame(l, columns=['Peptide', 'Genes', 'Affinity'])
df_pep_genes_aff.to_csv('../results/df_peptides_genes_aff_hla_a0201.tsv', index = False, sep = '\t')
