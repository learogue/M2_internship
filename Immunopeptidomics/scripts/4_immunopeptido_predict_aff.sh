#!/bin/bash
# ------------------------------------------------------------------------------------------------------------------
# Bash script : Prediction of peptide binding affinity to MHC class I molecules
# Author  : Léa ROGUE
# Date    : 28-04-2024
# Description : This script performs prediction of peptide binding affinity to MHC class I molecules using netMHC
# and selects strong binding peptides.
# ------------------------------------------------------------------------------------------------------------------

# Predict the binding affinity of the selected peptides to MHC class I molecules using NetMHCpan
alleles=(
    HLA-A0101
    HLA-A0201
    HLA-A0301
    HLA-B0702
    HLA-B0801
    HLA-B4402
    HLA-C0401
    HLA-C0501
    HLA-C0602
    HLA-C0701
    HLA-C0702
)
for allele in "${alleles[@]}"; do
    nohup netMHC -a "$allele" -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc/res_netmhc_${allele}_selected_blastp.out 2>&1 &
done

# Loop for each NetMHC output file
for file in ../results/res_netmhc/res_netmhc_*_selected_blastp.out; do
    # Extract file name
    base=$(basename "$file" .out)

    # Select th strong binding peptides
    grep 'SB' "$file" > "../results/netmhc_sb/${base}_sb.out"
done
