#!/bin/bash
# Search published shell-matrix proteomes against the Littorina proteome (run inside WSL).
set -e
B=~/mamba/envs/phylo/bin
D=/mnt/c/micah/shells/data
cd $D/shell_proteomes
for f in lottia_consensus crepidula_fornicata crepidula_atrasolea haliotis_laevigata pinctada_fucata cepaea_orfs; do
  $B/diamond blastp -q $f.faa -d $D/orthologs/lit -o $f.vs_lit.tsv --sensitive --max-target-seqs 5 --evalue 1e-10 --quiet \
    --outfmt 6 qseqid sseqid pident length evalue bitscore
done
wc -l *.vs_lit.tsv
