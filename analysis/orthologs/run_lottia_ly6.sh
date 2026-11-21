#!/bin/bash
B=~/mamba/envs/phylo/bin
cd /mnt/c/micah/shells/data/shell_proteomes
$B/diamond blastp -q lottia_ly6.faa -d ../orthologs/lit -o lottia_ly6.vs_lit.tsv --ultra-sensitive --max-target-seqs 10 --evalue 1 --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore
cat lottia_ly6.vs_lit.tsv
