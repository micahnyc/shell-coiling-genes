#!/bin/bash
# Search five annotated genomes directly for dispatched genes the annotation might have missed.
set -u
B=~/mamba/envs/phylo/bin
D=/mnt/c/micah/shells/data/absence_check
Q=/mnt/c/micah/shells/data/nautilus/disp_queries.faa
cd $D
for acc in GCF_028476545.2 GCF_000002075.1 GCF_034508935.2 GCF_046866305.1 GCF_016097555.1; do
  if [ ! -s $acc.fna ]; then
    curl -s -o $acc.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA"
    python3 -c "import zipfile,sys,shutil; z=zipfile.ZipFile(sys.argv[1]); n=[x for x in z.namelist() if x.endswith('_genomic.fna')][0]; shutil.copyfileobj(z.open(n), open(sys.argv[2],'wb'))" $acc.zip $acc.fna && rm $acc.zip
  fi
  [ -f $acc.mpi ] || $B/miniprot -t 8 -d $acc.mpi $acc.fna 2>/dev/null
  $B/miniprot -t 8 --gff --trans -N 5 $acc.mpi $Q > $acc.gff 2>/dev/null
  rm -f $acc.mpi
  echo "$acc done: $(grep -c mRNA $acc.gff) alignments"
done
