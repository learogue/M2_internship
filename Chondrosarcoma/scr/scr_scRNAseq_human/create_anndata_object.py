#!/usr/bin/env python3
# ----------------------------------------------------------------------------------------------------------------------------------------
# R script : Create raw data objects
# Auteur  : Léa ROGUE
# Date    : 31-03-2025
# Description : This script create anndata objects from 10X files.
# ----------------------------------------------------------------------------------------------------------------------------------------

# Import packages
import os
import scanpy as sc
import anndata as ad
from matplotlib import pyplot as plt

def create_anndata_object(source, data_type, output_dir):
    """
    Create AnnData object from 10X or Smart-seq2 data

    Input:
        source (str): The directory for 10X data
            10X : a directory containing 3 files
                barcodes.tsv
                    AAACGGGGTCTAGTGT 
                    AAAGCAATCGTTTATC
                    AACCGCGAGATATACG
                genes.tsv
                    ENSMUSG00000102693  4933401J01Rik
                    ENSMUSG00000064842  Gm26206
                    ENSMUSG00000102851  Gm18956
                matrix.tsv
                    %%MatrixMarket matrix coordinate integer general
                    55141 138 272219
                    415 1 1
        data_type (str): The type of the data :'10X'
        output_dir (str): The output dir
    
    Output:
        Anndata objects no filtered named object_number_of_dataset_samples.h5ad
            file.h5ad contains n_obs * n_vars = 90 * 55141 ; obs = rows = cells, var = columns = genes
                                obs: 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_MT', 'dataset'
                                var: 'MT', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'
   """
    # For 10X
    if data_type == '10X':
        # Create file names and paths
        name = source.replace('_10X', '')
        name_object = f'object_{name}_ori.h5ad'
        input_path = f'../../data/scrnaseq_data/{source}'

    # Create a directory for plots if it doesn't exist
    plot_dir = output_dir + f'/Plots/Plots_{name}'
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)

    if not os.path.isdir(output_dir + '/Objects/Objects_ori'):
        os.mkdir(output_dir + '/Objects/Objects_ori')

    # Reading the data
    if data_type == '10X':
        adata = sc.read_10x_mtx(input_path) # Read 10X data

    # Ensure unique var and obs names
    adata.var_names_make_unique() # Make gene names unique
    adata.obs_names_make_unique() # Make cell names unique

    # Calculate mitochondrial gene percentages
    adata.var['MT'] = adata.var_names.str.startswith('MT-') # Identify mitochondrial genes
    sc.pp.calculate_qc_metrics(adata, qc_vars = ['MT'], percent_top = None, log1p = False, inplace = True) # Calculate QC metrics

    # Add dataset name to the observation metadata
    adata.obs['dataset'] = adata.n_obs * [name]

    # Generate and save a violin plot for quality control metrics
    with plt.rc_context():
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_MT'], jitter = 0.4, size = 1.6, multi_panel = True, show = False)
        plt.savefig(f'{plot_dir}/violin_plot.pdf', bbox_inches = 'tight')

    # Save the AnnData object to a file
    adata.write(output_dir + f'/Objects/Objects_ori/{name_object}')

    print('Done')
