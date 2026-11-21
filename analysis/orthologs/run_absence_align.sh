#!/bin/bash
# Re-align with relaxed secondary-hit settings (-N 20 --outs=0.5), validated on Littorina LG17.
# Genomes are downloaded by run_absence_check.sh.
set -u
B=~/mamba/envs/phylo/bin
D=/mnt/c/micah/shells/data/absence_check
Q=/mnt/c/micah/shells/data/nautilus/disp_queries.faa
cd $D
for acc in GCF_028476545.2 GCF_000002075.1 GCF_034508935.2 GCF_046866305.1 GCF_016097555.1; do
  $B/miniprot -t 8 --gff --trans -N 20 --outs=0.5 $acc.fna $Q > $acc.gff 2>/dev/null
  echo "$acc: $(grep -c mRNA $acc.gff) alignments"
done
