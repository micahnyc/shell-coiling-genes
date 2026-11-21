#!/bin/bash
# Align dispatched and diaphanous proteins to the Nautilus pompilius genome (GCA_018389105.1).
set -e
B=~/mamba/envs/phylo/bin
cd /mnt/c/micah/shells/data/nautilus
[ -f nautilus.mpi ] || $B/miniprot -t 8 -d nautilus.mpi nautilus.fna
$B/miniprot -t 8 --gff --trans -N 3 nautilus.mpi queries.faa > nautilus_hits.gff
grep -c "mRNA" nautilus_hits.gff
