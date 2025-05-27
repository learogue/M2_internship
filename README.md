# IMoST M2 Bioinformatic internship

This repository contains scripts created during my M2 internship in Bioinformatics, focused on the identification of therapeutic targets for radiopharmaceutical therapy in chondrosarcoma.

The objective of this work is to characterize tumor-specific markers in chondrosarcoma. The project is based on transcriptomic and immunopeptidomic data analyses to identify potential peptides present on the surface of cancerous chondrosarcoma cells. We focused on genes classified as Cancer-Testis Antigens (CTA), which are expressed in many cancers. Specifically, we explored the expression of these CTA genes in chondrosarcoma and their potential to be presented on the surface of tumor cells via the HLA-I complex, through an immunopeptidomic analysis.

The first part concern transcriptomic analysis with public dataset of chondrosarcoma (www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7264), healthy tissues (gtexportal.org) dataset which is single cell RNA-seq (GSE184118).
The 2nd part concer immunopeptidomic analysis with the Immune Epitope Database exploration (www.iedb.org).


## Usage and requirements
This project was developed and executed on a virtual machine provided by IFB Biosphère (8 CPUs, 32 GB RAM) and on local. It uses multiple languages including R, Python, Bash, and MySQL. All required packages are listed in the cta.yaml file and can be installed using Conda.

We also used NetMHC 4.0 (services.healthtech.dtu.dk/services/NetMHC-4.0) for HLA-I binding predictions (authorization is required to use this tool).

## Installation
```
git clone https://github.com/learogue/M2_internship.git
cd M2_internship 
conda env create -f env/cta.yaml
```
You can execute the scripts directly or use RStudio to run the RMarkdown files.

## Contents
```
├── Chondrosarcoma
│   ├── data
│   │   ├── CTA_list_clean.txt                          # Complete list of CTA (including validated and putatives)
│   │   ├── MHC_genes.txt                               # MHC-I and MHC-II genes coding
│   │   ├── immune_cells_genes.tsv                      # Signature genes per immune cells type (Bindea et al., 2013)
│   │   └── list_files_used_publi.txt                   # Files list because there is duplicates unassociated patient data
│   ├── results
│   └── scripts
│       ├── 00_chondro_download_data.R                  # Download microarray data (chondrosarcoma)
│       ├── 01_clean_files.sh                           # Remove duplicates files
│       ├── 02_chondro_preprocess.R                     # Create raw and normalized (RMA) objects
│       ├── 03_chondro_qc.R                             # QC
│       ├── 04_chondro_processing.R                     # Matrixes constructions (genes, immune cells)
│       ├── 05_chondro_expression_analysis.Rmd          # Transcriptomic analysis
│       ├── 06_chondro_survival_analysis.Rmd            # Survival analysis
│       ├── 07_expression_cta_normal_tissues_gtex.R     # Download data and matrix construction from GTEx to have mean TPM per CTA
│       ├── 08_normal_tissues_cta_expression.Rmd        # Explore CTA expression in normal tissues and integrate CTA expression in Chondrosarcoma
│       ├── 09_pseudo_bulk_compute_tpm.py
│       ├── 10_pseudo_bulk_scrnaseq.Rmd
│       ├── 11_sc_analysis_conv_chondro.ipynb
│       ├── pdf_html_md                                 # Dir containing pdf and html files with figures from Rmarkdown
│       ├── scr_scRNAseq_hum                            # Contain scripts to pre-process scRNAseq data
│       │   ├── README.md                               # Explain how to use these scripts
│       │   ├── apply_filters.py                        # Filter genes and cells and QC
│       │   ├── create_anndata_object.py                # Create objects
│       │   ├── main.py                                 # Main script to execute the others
│       │   └── merge.py                                # Merge differents objects
│       └── tex_files                                   # Contain a file to genereate list of figures in pdf (Rmarkdown)
├── Immunopeptidomics
│   ├── data
│   ├── results
│   └── scripts
│       ├── 0_immunopeptido_extract_data_iedb.sql       # Extract data from SQL database
│       ├── 1_immunopeptido_create_peptides_table.py    # Peptide tables from SQL table (difficulties to explore database with mySQL)
│       ├── 2_immunopeptido_processing.sh               # Peptides BLASTp on selected genes sequences 
│       ├── 3_immunopeptido_fasta_aligned_pep.py        # Create fasta of aligned peptides and select hits with 0 mismatchs
│       ├── 4_immunopeptido_predict_aff.sh              # Predict affinity with netMHC tool
│       ├── 5_immunopeptido_analysis.py                 # Explore predict data
│       └── 6_immunopeptido_plot.Rmd                    # Integrate predicted data with expression
├── LICENSE
├── README.md
└── cta.yaml                                            # yaml file to generate conda environment
```

## Licence
Creative Commons Legal Code - CC0 1.0 Universal