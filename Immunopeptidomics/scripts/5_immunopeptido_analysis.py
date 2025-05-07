#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------
# Python script : Analyze the number of peptides per gene
# Author  : Léa ROGUE
# Date    : 02-05-2025
# Description : This script analyzes the number of peptides per gene from the results of the strong binding peptides for
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

# Count the number of peptides per gene and adding column with genes in the file with genes
d_nb_pep_gene = defaultdict(int)
with open('../results/netmhc_sb/res_netmhc_HLA-A0201_selected_blastp_sb_annotated.out', 'w') as f_w:
    with open('../results/netmhc_sb/res_netmhc_HLA-A0201_selected_blastp_sb.out', 'r') as f_r:
        for lig in f_r:
            genes = str(d_pep_gene[lig.strip().split()[10]]).replace('\'', '').replace('[', '').replace(']', '')
            d_nb_pep_gene[genes] += 1
            lig = lig.strip() + '\t' + genes
            f_w.write(lig + '\n')

# Sorted dictionary
sorted_dict = dict(sorted(d_nb_pep_gene.items(), key = lambda item: item[1], reverse = True))

# Create dataframe
df = pd.DataFrame(list(sorted_dict.items()), columns = ['Genes', 'Peptides number'])

# Save
df.to_csv('../results/nb_pep_per_gene_sb_hla_a0201.tsv', index = False, sep = '\t')

# Count the number of peptides per unique and take affinity
d_pep_id_affinity = defaultdict(list)
d_nb_pep_gene = defaultdict(int)
with open('../results/netmhc_sb/res_netmhc_HLA-A0201_selected_blastp_sb.out', 'r') as f:
    for lig in f:
        lig = lig.strip().split()
        for gene in d_pep_gene[lig[10]]:
            d_nb_pep_gene[gene] += 1
        d_pep_id_affinity[lig[10]].append(float(lig[12]))

# Gene length
d_genes_length = defaultdict(int)
in_seq = False
with open('../data/proteine_seq_targeted_cta.fasta') as f:
    for lig in f:
        lig = lig.strip()
        
        if lig.startswith('>'):
            in_seq = False
            gene_name = re.search(r'GN=([\w-]+)', lig)
            gene_name = gene_name.group(1)
            d_genes_length[gene_name] = ''
        else:
            in_seq = True

        if in_seq:
            d_genes_length[gene_name] += lig.strip()

# Sorted dict
sorted_dict = dict(sorted(d_nb_pep_gene.items(), key = lambda item: item[1], reverse = True))

# Create dataframe
df = pd.DataFrame(list(sorted_dict.items()), columns = ['Genes', 'Nb_peptides_per_gene'])

# Calculate the length of each gene
for gene in d_genes_length:
    d_genes_length[gene] = len(d_genes_length[gene])

# Adding protein length to the dataframe
df['Protein_length'] = df['Genes'].map(d_genes_length)

# Scatter plot with protein length and number of peptides
plt.figure(figsize = (10, 6))
sns.regplot(
    data = df,
    x = 'Protein_length',
    y = 'Nb_peptides_per_gene',
    scatter = True,
    line_kws = {'color': 'red'},
    ci = None)  # Delete confidence interval
plt.title('Number of peptides per gene vs protein length')
plt.xlabel('Protein length')
plt.ylabel('Number of peptides per gene')
plt.savefig('../results/nb_pep_per_gene_sb_hla_a0201.png')
plt.show()

# Normalize the number of peptides per gene by the protein length
df['Pep_normalized'] = df['Nb_peptides_per_gene'] / df['Protein_length']

# Read the data
df_int = pd.read_csv('../../Chondrosarcoma/results/whole_gene_int_CTA_sign_imm_clean.tsv', sep = '\t')

# Select CTA genes
df_int = df_int[df_int['SYMBOL'].isin(df['Genes'])]

# Mean expression level
df_int['Mean_expression'] = df_int.iloc[:,3:].mean(axis = 1)

# Select data for the scatter plot
df_int = df_int[['SYMBOL', 'Mean_expression']].copy()
df_int = df_int.merge(df, left_on = 'SYMBOL', right_on = 'Genes', how = 'inner')
df_int = df_int[['SYMBOL', 'Mean_expression', 'Pep_normalized']].copy()

# Scatter plot with number of peptides and expression level in chondrosarcoma
plt.figure(figsize = (10, 6))
sns.scatterplot(data = df_int, x = 'Mean_expression', y = 'Pep_normalized')
plt.title('Normalized umber of peptides per gene vs expression level in chondrosarcoma')
plt.xlabel('Expression level in chondrosarcoma')
plt.ylabel('Normalized number of peptides per gene')
plt.savefig('../results/scatter_expr_nb_pep_norm_hla_a0201.png')
plt.show()

# Prepare data for scatter plot with expression level in the tumor and affinity
d_gene_affinity = defaultdict(list)
for pep, l_aff in d_pep_id_affinity.items():
    for genes in d_pep_gene[pep]:
        d_gene_affinity[genes] += l_aff

# Take the minimum affinity
for gene in d_gene_affinity:
    d_gene_affinity[gene] = min(d_gene_affinity[gene])

# Create dataframe
df_affinity = pd.DataFrame(list(d_gene_affinity.items()), columns = ['Genes', 'Affinity'])

# Merge the dataframes
df_int = df_int.merge(df_affinity, left_on = 'SYMBOL', right_on = 'Genes', how = 'inner')
df_int = df_int[['SYMBOL', 'Mean_expression', 'Pep_normalized', 'Affinity']].copy()

# Save df with affinity
#df_int[['SYMBOL', 'Affinity']].to_csv('../results/df_min_affinity_expr_hla_a0201.tsv', index = False, sep = '\t')

# Create a scatter plot with mean affinity and expression level in chondrosarcoma
plt.figure(figsize=(10, 6))
sns.scatterplot(data = df_int, x = 'Affinity', y = 'Mean_expression')
plt.xscale('log')
plt.title('Minimum affinity per gene vs expression level in chondrosarcoma')
plt.xlabel('Minimum affinity (nM)')
plt.ylabel('Expression level in chondrosarcoma')
plt.savefig('../results/scatter_affinity_expr_hla_a0201.png')

# Create interactive figure
fig = go.Figure(data=go.Scatter(
    x=df_int['Affinity'],
    y=df_int['Mean_expression'],
    mode='markers',
    text=df_int['SYMBOL']))

# Add title and labels
fig.update_layout(
    title='Minimum affinity per gene vs expression level in chondrosarcoma',
    xaxis_title='Affinity (nM)',
    yaxis_title='Expression level in chondrosarcoma',
    hovermode='closest',
    xaxis=dict(type='log')) # Log scale for x-axis

# Save the figure
fig.write_html('../results/scatter_affinity_expr_hla_a0201_interactive.html')



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

