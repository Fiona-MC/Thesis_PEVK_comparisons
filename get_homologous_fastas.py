# Author: Fiona
# Date: 9 March 2021

# this is meant to be directly downstream of the reciprocal_blast_exonL.R file -- matched exons is the csv file
# this code takes a table of coordinates in 12345:12456 form in a csv file and outputs fasta files with
# each homologous exon -- then puts them in a folder that is called homologous_fasta
# expects a file as input that has columns 2-n with titles that are species and that there are .fasta files with the titin sequence from each species called ./Species_name_ttn.fasta

# to run: python3 get_homologous_fastas.py -i matched_exon_data.csv 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein

import argparse
import os
import csv


def read_csv_file(csv_filename):
    ''' Takes the filename of csv file and reads it into a table with each row being a single homologous exon
    Returns the coordinates table
    '''
    exon_coordinate_table = []
    with open(csv_filename, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        for row in reader:
            # check if we've reached the end of the real rows (ie this row is completely filled with NA's)
            if row[2:] == ["NA" for i in range(len(row)-2)]:
                break
            # If it has real content, put into the table
            exon_coordinate_table.append(row)
    return exon_coordinate_table


def makeFasta(exonCoordL, speciesL, inFolder, outpath = "./"):
    ''' Writes fasta from a list of coordinates
    Takes a single row of the exon_coordinate_table and the list of species (from the first row of the exon coord table)
    writes a fasta file with the title being the content of the first entry
    (name of the human exon)
    and the subsequent entries being the sequences from the various species that correspond to this exon
    '''
    # initialize list of exons (this will contain seqRecord objects)
    allExons = []

    # Use coordinates list to get all of the sequences from the ttn sequences
    for i in range(1,len(exonCoordL)):
        exonCoords = exonCoordL[i]
        # If this has an entry
        if exonCoords != "NA":
            # get the fasta sequence of the whole titin molecule 
            thisSpecies = speciesL[i]
            thisSpeciesFastaName = inFolder+thisSpecies+"_ttn.fasta"
            seq_record = SeqIO.read(thisSpeciesFastaName, "fasta") 
            exonStart = int(exonCoords.split(":")[0])
            exonEnd = int(exonCoords.split(":")[1])

            thisSpeciesExon = Seq(str(seq_record.seq[exonStart:exonEnd+1]), generic_dna)
            nt_seq_record = SeqRecord(thisSpeciesExon, id=thisSpecies+ "_" + str(exonStart) + ":" + str(exonEnd), description='')
            allExons.append(nt_seq_record)

    # write the new fasta with all of the homologs found
    newFastaName = outpath+"exon_"+exonCoordL[0]+"_homologs.fasta"
    SeqIO.write(allExons, newFastaName, "fasta")
    return 


def makeAllFastas(exon_coordinate_table, inFolder, outpath):
    ''' Runs through all homologous exons
    and writes fasta files for them
    '''
    speciesL = exon_coordinate_table[0]
    for exonCoordL in exon_coordinate_table:
        if exonCoordL != speciesL: # If this isn't the list of species (ie first row)
            # Write fasta file for this exon
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            makeFasta(exonCoordL, speciesL, inFolder, outpath = outpath)
    return


def main():
    # get filename
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                   help="file name of CSV",
                   type=str, default = "matched_exon_data.csv")
    ap.add_argument("-d", "--inFolder",
                   help="file folder name for _ttn.fasta files",
                   type=str, default = "./")

    args = ap.parse_args()
    fileName = args.input
    inFolder = args.inFolder

    exon_coordinate_table = read_csv_file(fileName)
    makeAllFastas(exon_coordinate_table, inFolder, "./homologous_fasta/")
    

if __name__=="__main__":
    main()