#!/bin/bash
# Align dispatched and diaphanous proteins to the Nautilus pompilius genome (GCA_018389105.1).
set -e
B=~/mamba/envs/phylo/bin
cd /mnt/c/micah/shells/data/nautilus
[ -f nautilus.mpi ] || $B/miniprot -t 8 -d nautilus.mpi nautilus.fna
$B/miniprot -t 8 --gff --trans -N 20 --outs=0.5 nautilus.mpi disp_queries.faa > ../absence_check/nautilus.gff
grep -c "mRNA" nautilus_hits.gff
