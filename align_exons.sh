#!/bin/bash
# Author: Fiona Callahan
# This assumes you are on the server where clustalo exists

cd ./homologous_fasta/

for infileName in ./*
do
outname=$(echo "${infileName%.*}_aligned.fasta")

/programs/clustalo --infile $infileName --threads 8 --MAC-RAM 8000 --verbose  --outfmt fa --outfile $outname --output-order tree-order --seqtype dna
done

mkdir ../aligned_exons
mv *_aligned.fasta ../aligned_exons/

cd ..
