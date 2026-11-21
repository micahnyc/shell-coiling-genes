#!/bin/bash
# Align, trim and build ML trees for each family (run inside WSL Ubuntu).
# Tools: micromamba env ~/mamba/envs/phylo (mafft 7.526, trimal, iqtree3).
set -e
B=~/mamba/envs/phylo/bin
cd /mnt/c/micah/shells/data/gene_trees
for fam in formin dispatched; do
  $B/mafft --auto --thread 4 $fam.faa > $fam.aln.faa 2> $fam.mafft.log
  $B/trimal -in $fam.aln.faa -out $fam.trim.faa -automated1
  $B/iqtree3 -s $fam.trim.faa -m MFP -mset LG,WAG,JTT -B 1000 -T AUTO --threads-max 6 -redo --prefix $fam > $fam.iqtree.log 2>&1
  grep -E "Best-fit model|Log-likelihood of consensus" $fam.iqtree || true
done
