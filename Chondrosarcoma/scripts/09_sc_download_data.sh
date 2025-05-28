#!/bin/bash
# -----------------------------------------------------------------------------------------------------------------
# bash script : Download processed scRNAseq  files and rename
# Auteur  : Léa ROGUE
# Date    : 28-03-2025
# Description : This script downloads processed scRNA-seq data files, extracts them from compressed archives, 
# removes the original compressed files, and renames the extracted directories for easier downstream analysis.
# -----------------------------------------------------------------------------------------------------------------

# Change dir
cd ../data/scrnaseq_data/

# Dowload the data
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184118/suppl/GSE184118_RAW.tar

# Untar
tar -xvf GSE184118_RAW.tar

# Untar dir
for file in GSM557818*.tar.gz; do
  tar -xzvf "$file"
done

# Rmove tar files
rm -rf GS*

# Old dir names
DIR_NAMES=("L07" "L28" "L31" "L43" "L44" "L49" "L63" "L75" "L80" "L81" "L83")

# New dir names
NEW_NAMES=(
  "1_Low_L07_10X"
  "2_Low_L28_10X"
  "3_High_L31_10X"
  "4_Cos_L43_10X"
  "5_High_L44_10X"
  "6_Ben_L49_10X"
  "7_Med_L63_10X"
  "8_FF_L75_10X"
  "9_Med_L80_10X"
  "10_Low_L81_10X"
  "11_Ded_L83_10X"
)

# Loop to rename
for i in "${!DIR_NAMES[@]}"; do
	mv "${DIR_NAMES[$i]}" "${NEW_NAMES[$i]}"
done

