#bin/bash/env python3
# 15/04/2025

from collections import defaultdict
import pandas as pd

# Take disease id and disease name
d_disease = defaultdict()
with open('../results/db_tables/disease.tsv', 'r') as f:
    for lig in f:
        lig = lig.strip().split('\t')
        d_disease[lig[0]] = lig[1]

# take tcell line with a cancer disease
d_tcell = defaultdict(list)
with open('../results/db_tables/t_cell.tsv', 'r') as f:
    for lig in f:
        lig = lig.strip().split('\t')

        disease_name = []
        for col in [3, 4, 5]:
            key = lig[col]
            if key in d_disease:
                label = d_disease[key]
                if label not in disease_name:
                    disease_name.append(label)

        if disease_name:
            d_tcell[lig[0]] = [lig[0], lig[1], lig[2]] + [", ".join(disease_name), "t_cell"]

# take bcell line with a cancer disease
d_bcell = defaultdict(list)
with open('../results/db_tables/b_cell.tsv', 'r') as f:
    for lig in f:
        lig = lig.strip().split('\t')

        disease_name = []
        for col in [2, 3, 4]:
            key = lig[col]
            if key in d_disease:
                label = d_disease[key]
                if label not in disease_name:
                    disease_name.append(label)

        if disease_name:
            d_bcell[lig[0]] = [lig[0], lig[1]] + [", ".join(disease_name), "b_cell"]

# MHC epitope
d_mhc = defaultdict(list)
with open('../results/db_tables/mhc_epitope.tsv', 'r') as f:
    for lig in f:
        lig = lig.strip().split('\t')
        disease_labels = []

        if lig[1] in d_disease.keys():
            d_mhc[lig[0]] = [lig[0], lig[2], lig[3]] + [d_disease[lig[1]], "mhc_epitope"]

# Create dataframe
df_bcell = pd.DataFrame.from_dict(d_bcell, orient='index', columns=['id', 'mhc_type', 'disease', 'source'])
df_mhc = pd.DataFrame.from_dict(d_mhc, orient='index', columns=['id', 'mhc_type', 'mhc_allele', 'disease', 'source'])
df_tcell = pd.DataFrame.from_dict(d_tcell, orient='index', columns=['id', 'mhc_type', 'mhc_allele', 'disease', 'source'])

# Merge
df_all = pd.concat([df_bcell[['id', 'mhc_type', 'disease', 'source']].assign(source='b_cell'),
                    df_mhc[['id', 'mhc_type', 'mhc_allele', 'disease', 'source']].assign(source='mhc_epitope'),
                    df_tcell[['id', 'mhc_type', 'mhc_allele', 'disease', 'source']].assign(source='t_cell')])

# Merge same ids
df_merge = df_all.groupby('id').agg(lambda x: ', '.join(str(i) for i in x.unique())).reset_index()

# Take sequences in object_seq thanks to id_object link to epitope id in epitope_object.tsv
d_obj_seq = defaultdict()
d_epi_seq = defaultdict()
with open('../results/db_tables/epitope_object.tsv', 'r') as f1:
    with open('../results/db_tables/object_seq.tsv', 'r') as f2:
        # Take ids with a sequences
        for lig in f2:
            lig = lig.strip().split('\t')

            if lig[1] != '\\N':
                d_obj_seq[lig[0]] = lig[1]

        # Take epitope ids with a sequence
        for lig in f1:
            lig = lig.strip().split('\t')

            if lig[1] in d_obj_seq.keys():
                d_epi_seq[lig[0]] = d_obj_seq[lig[1]]

# Create df and transform id col as character
df_epi_seq = pd.DataFrame(list(d_epi_seq.items()), columns = ['id', 'sequence'])
df_epi_seq['id'] = df_epi_seq['id'].astype(str)
df_merge['id'] = df_merge['id'].astype(str)

# Merge sequences with ids
df_merge = df_merge.merge(df_epi_seq, on = 'id', how = 'inner')

# Save
df_merge.to_csv("../results/table_peptides_cancer_human.tsv", sep="\t", index=False)
