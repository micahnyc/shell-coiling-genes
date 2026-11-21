#!/bin/bash
set -e
B=~/mamba/envs/phylo/bin
cd /mnt/c/micah/shells/data/dreissena_sc
[ -f dre.dmnd ] || $B/diamond makedb --in dreissena_proteins.faa -d dre --quiet
$B/diamond blastp -q ${1:-queries.faa} -d dre -o ${2:-hits.tsv} --sensitive --max-target-seqs 20 --evalue 1e-10 --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore
