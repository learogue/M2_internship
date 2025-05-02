#!/bin/bash
# ------------------------------------------------------------------------------------------------------------------
# Bash script : Prediction of peptide binding affinity to MHC class I molecules
# Author  : Léa ROGUE
# Date    : 28-04-2024
# Description : This script performs prediction of peptide binding affinity to MHC class I molecules using netMHC
# and selects strong binding peptides.
# ------------------------------------------------------------------------------------------------------------------

# Predict the binding affinity of the selected peptides to MHC class I molecules using NetMHCpan
nohup netMHC -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_selected_blastp.out 2>&1 &

# Select th strong binding peptides
cat res_netmhc_selected_blastp.out | grep 'SB' > strong_binding_peptides_hla_a0201.out

# Other MHC class I alleles
nohup netMHC -a HLA-A0101 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_a0101_selected_blastp.out 2>&1 &
nohup netMHC -a HLA-A0301 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_a0301_selected_blastp.out 2>&1 &
nohup netMHC -a HLA-B0702 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_b0702_selected_blastp.out 2>&1 &
nohup netMHC -a HLA-B0801 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_b0801_selected_blastp.out 2>&1 &
nohup netMHC -a HLA-B4402 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_b4402_selected_blastp.out 2>&1 &


nohup netMHC -a HLA-C0401 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_c0401_selected_blastp.out 2>&1 &
nohup netMHC -a HLA-C0501 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_c0501_selected_blastp.out 2>&1 &
nohup netMHC -a HLA-C0602 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_c0602_selected_blastp.out 2>&1 &
nohup netMHC -a HLA-C0701 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_c0701_selected_blastp.out 2>&1 &
nohup netMHC -a HLA-C0702 -l 8,9,10 ../results/seq_pep_aligned.fasta > ../results/res_netmhc_hla_c0702_selected_blastp.out 2>&1 &