# 15-04-2024

# Make fasta file with peptide sequences
awk -F'\t' 'NR > 1 {print ">" $1 "\n" $6}' results/table_peptides_cancer_human.tsv > results/seq_pep.fasta

# Make db to use blast with selected genes
makeblastdb -in data/proteine_seq_targeted_cta.fasta -dbtype prot -out data/db/targeted_cta_db

# Align petides with sequences of selected genes
blastp -query results/seq_pep.fasta -db data/db/targeted_cta_db -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length bitscore evalue pident mismatch' -evalue 0.01 -num_threads 6 -out results/results_blastp_mismatch_eval_0_01.tsv

# Select hits with 0 or 1 mismatchs
awk -F'\t' '($13 == 0 || $13 == 1)' results/results_blastp_mismatch_eval_0_01.tsv > results/selected_results_blastp_mismatch_eval_0_01.tsv
