# Thesis_PEVK_comparisons

Instructions to run full pipeline for thesis. Steps 1-2 find PEVK exons. Steps 3-6 find homologous exons. Steps 7-10 run a likelihood ratio test on different evolutionary models using PAML.
This code depends on a lot of different packages -- notably: PEVK_finder depends on python2 and the rest of the python files are written in python3.

All of these instructions assume you are in the folder with these instructions and all of the code, and that the outputs of previous steps have not been moved out of that folder.


1. get the ttn fasta sequences from NCBI and put them in a folder with files called Species_name_ttn.fasta 
Put all of them in a folder called Data
Note: Make sure the species name on the file exactly matches the one in the FASTA header.

python3 find_gaps.py -i /Path/input/file/Species_name_ttn.fasta
will find the gap length and coordinates of gaps

2. Run Pevk_finder_v_1.py on all files 
python2 -W ignore ./pevk_finder_public-master/pevk_finder_v_1.py -i /input/folder/with/sequences/ -w 10 -r 0.54 -l 12 -p True -o /output/folder

Note: The following steps assume the output folder was the folder with all of this code in it.

From these, you can count the exons and total exon length using get_pevk_stats.py

python3 get_pevk_stats.py -i ./exon_lengths_and_ratios/ -l testLabel -o ./

3. Make blast DBs (make sure to add the location of the databases to your path)
Run the following at the command line:
# Make Databases
for exonfile in *.fasta
do
# Extract query species name
species_name=${exonfile%_PEVK_exons_NT_unbounded_10_0.54_12.fasta}
makeblastdb -in $exonfile -out /Path/to/output/blast/db/$species_name"_db" -dbtype nucl
done

4. Make reciprocal blast csv files
run the code in reciprocal_blast.sh
Note: If you want to run it from the command line:
chmod u+x reciprocal_blast.sh
./reciprocal_blast.sh

Note: careful about running this more than once --  it will duplicate and overwrite things unexpectedly -- if you need to run it again, remove the previous outputs first.

5. Match the exons in CSV format
Run the code in reciprocal_blast_exonL.R 
(will need a few edits for the folder directories in the file)
The output file is called /nt_matches_csv/matched_exon_data.csv
Move it out of that folder

6. Get actual exons and put them in folder of fastas
Run the code in get_homologous_fastas.py
ie at command line
python3 get_homologous_fastas.py -i matched_exon_data.csv -d ./Data/

7. align all exons: -- do this on the server because it has clustal omega installed
chmod u+x ./align_exons.sh
./align_exons.sh

8. make tree for all species using timetree.org 
You can use the file called SpeciesL.txt as input for that file
Note: a few of the species have multiple possible names (eg Neomonachus vs Monachus) and that messes up subsequent steps so be careful of that.

9. If you want to run the PAML codeml LRT, you need tree, phylip file with seqs, and codeml.ctl control file. There are samples of these in the PAML_sample_inputs folder.


Make phylip file:
python3 get_phylip_all_exons.py

FOR phylip file--to get it prepped for PAML
1) add one extra space (for a total of 2 spaces) between each species name and the sequence
2) add I for "interlaved" after the two numbers on the first line

for TREE to get prepped for PAML
1) change all names to have only the first 10 chars
2) make sure the tree has a semicolon at the end
3) add branch labels ie $1 or #1 to marine/subterr clades -- see PAML manual

ctl file: for H0, clock=0, model=0; for H1 clock=0, model=2

10. Run LRT for H2 vs H0 and H1 vs H0 -- calculate test stat (by hand, 2*(loglikelihoodH1-loglikelihoodH0)) and run 

testStat=
chi2 1 $testStat
