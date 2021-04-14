#!/bin/bash
# Fiona adapted this from Kathleen's code
# this code blasts every species exons against every other species in the folder you're in
# assumes you have already made the exon databases
# assumes you are in the folder with the original fasta files, which then contains output of pevk finder (ie ./Data/unbounded/untranslated)

# for full pipeline see readme_reciprocal_blast.txt

# get species list:
rm speciesL.txt  # this line makes sure that we aren't just adding more species to list
touch speciesL.txt
cd ./Data/
for f in *.fasta
do
    species=${f%_ttn.fasta}
    echo $species >> ../speciesL.txt
done
cd ..

# Loop through all fasta files containing the exon seqs
for f in ./unbounded/untranslated/*.fasta
do
    # Extract query species name
    species_name=${f%_PEVK_exons_NT_unbounded_10_0.54_12.fasta}

    # Loop through all species and blast against db created out of that specis's exons
    while read -r line
    do
	s="$line"
	blastn -evalue 1e-6 -query $f -db $s"_db" -out $species_name"_and_"$s"_blast" -outfmt 6
	mv $species_name"_and_"$s"_blast" .
    done < speciesL.txt
done

# Translate all BLAST inputs to CSV files so they can be read by R
for f in *blast
do
    # Extract blast identifier
    blast_identifier=${f%_blast}

    # Convert file in csv
    tr '\t' ',' < $f > $blast_identifier"_RB.csv"
done

mkdir nt_matches

# Move blast files into separate folder
mv *blast nt_matches

mkdir nt_matches_csv

# Move csv files into separate folder
mv *_RB.csv nt_matches_csv


