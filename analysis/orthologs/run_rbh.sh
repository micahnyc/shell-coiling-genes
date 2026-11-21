#!/bin/bash
# Reciprocal best hits of Littorina query proteins in three species (run inside WSL).
set -e
B=~/mamba/envs/phylo/bin
D=/mnt/c/micah/shells/data
cd $D/orthologs
OPT="--sensitive --max-target-seqs 5 --evalue 1e-10 --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore"
[ -f lit.dmnd ] || $B/diamond makedb --in $D/littorina_v2/protein.faa -d lit --quiet
[ -f myt.dmnd ] || $B/diamond makedb --in GCF_965363235.1_xbMytGall1.hap1.1_protein.faa.gz -d myt --quiet
[ -f cvi.dmnd ] || $B/diamond makedb --in GCF_053477285.1_ASM5347728v1_protein.faa.gz -d cvi --quiet
cp -n $D/dreissena_sc/dre.dmnd . 2>/dev/null || true
for sp in myt cvi dre; do
  $B/diamond blastp -q littorina_queries.faa -d $sp -o fwd_$sp.tsv $OPT
  cut -f2 fwd_$sp.tsv | sort -u > hits_$sp.txt
  # extract hit proteins from the target fasta, then search back against all Littorina proteins
  case $sp in
    myt) src="zcat GCF_965363235.1_xbMytGall1.hap1.1_protein.faa.gz";;
    cvi) src="zcat GCF_053477285.1_ASM5347728v1_protein.faa.gz";;
    dre) src="cat $D/dreissena_sc/dreissena_proteins.faa";;
  esac
  $src | tr -d '\r' | awk 'NR==FNR{w[$1];next} /^>/{split(substr($0,2),a," "); k=(a[1] in w)} k' hits_$sp.txt - > hits_$sp.faa
$/,"")} NR==FNR{w[$1];next} /^>/{split(substr($0,2),a," "); k=(a[1] in w)} k' hits_$sp.txt - > hits_$sp.faa
  $B/diamond blastp -q hits_$sp.faa -d lit -o rev_$sp.tsv $OPT
done
wc -l fwd_*.tsv rev_*.tsv
