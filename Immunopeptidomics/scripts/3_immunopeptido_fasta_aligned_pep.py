#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------
# Python script : Create a fasta file with aligned peptides
# Author  : Léa ROGUE
# Date    : 21-04-2025
# Description : This script permit to create a fasta files with peptide sequences that are aligned by blastp to run the
# netMHC tool to predict affinity of peptides to MHC.
# ----------------------------------------------------------------------------------------------------------------------

# Import packages
from collections import defaultdict

# Read gene names and their corresponding entry names
d_genes = defaultdict()
with open('../data/entry_name_gene_name.tsv', 'r') as f:
    for lig in f:
        lig = lig.strip().split()
        d_genes[lig[0]] = lig[1]

# Read peptides IDs and their corresponding gene
d_pep_gene_id = defaultdict(list)
pep_id = set()
with open('../results/selected_results_blastp_0_mismatch.tsv') as f:
    for lig in f:
        pep_id.add(lig.strip().split()[0])

# Extract peptides from the fasta file
d_pep_aligned = defaultdict(str)
in_seq = False
with open('../results/seq_pep.fasta') as f:
    for lig in f:
        lig = lig.strip()
        
        if lig.startswith('>'): # if line starts with '>', verify if it is an aligned peptide
            in_seq = False
            if lig.replace('>', '') in pep_id:# if it is an aligned peptide, save the id and start saving the sequence
                id = lig
                d_pep_aligned[id] = ''
                in_seq = True
        elif in_seq:
            d_pep_aligned[id] += lig

# Save the aligned peptides in a new fasta file
with open('../results/seq_pep_aligned.fasta', 'w') as f:
    for id, seq in d_pep_aligned.items():
        f.write(id + '\n' + seq + '\n')
