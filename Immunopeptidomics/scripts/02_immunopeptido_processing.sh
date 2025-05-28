#!/bin/bash
# ------------------------------------------------------------------------------------------------------------------
# Bash script : Peptide sequence alignment with selected CTA protein sequences using BLAST
# Author  : Léa ROGUE
# Date    : 15-04-2024
# Description : This script performs several steps to align cancer-associated human peptides with a set of selected 
# cancer-testis antigen (CTA) protein sequences. First, it extracts peptide sequences from TSV file and converts 
# them into a FASTA format. A BLAST database is then created from the CTA protein sequences. The peptide sequences 
# are aligned against this database using blastp. The results are filtered to retain only those alignments with 0
# or 1 mismatches. Finally, gene names are extracted from the CTA fasta file to allow mapping of hits to gene 
# identifiers.
# ------------------------------------------------------------------------------------------------------------------

# Make fasta file with peptide sequences
awk -F'\t' 'NR > 1 {print ">" $1 "\n" $6}' results/table_peptides_cancer_human.tsv > results/seq_pep.fasta

# Make db to use blast with selected genes
makeblastdb -in data/proteine_seq_targeted_cta.fasta -dbtype prot -out data/db/targeted_cta_db

# Align petides with sequences of selected genes
blastp -query results/seq_pep.fasta -db data/db/targeted_cta_db -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length bitscore evalue pident mismatch' -num_threads 6 -out results/results_blastp.tsv

# Select hits with 0 mismatch
awk -F'\t' '($13 == 0)' results/results_blastp.tsv > results/selected_results_blastp_0_mismatch.tsv

# Take gene names and sseid
grep '^>' proteine_seq_targeted_cta.fasta | tr -d '>' | awk '{
    for (i=1; i<=NF; i++) {
        if ($i ~ /^GN=/) {
            sub(/^GN=/, "", $i);
            print $1, "\t", $i;
        }
    }
}' > entry_name_gene_name.tsv 
