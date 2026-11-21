#!/bin/bash
# Positive control: Littorina LG17 carries 6 lineage genes and DISP1 in its annotation.
B=~/mamba/envs/phylo/bin
D=/mnt/c/micah/shells/data
$B/miniprot -t 4 --gff --trans -N 20 --outs=0.5 $D/littorina_v2/NC_090261.1.fa $D/nautilus/disp_queries.faa > $D/absence_check/control_LG17.gff 2>/dev/null
grep -c mRNA $D/absence_check/control_LG17.gff
