# Author: Fiona Callahan
# this file goes through a folder of homologous fastas that are pre-aligned -- see PAML NOTES word doc, 
# checks that they have a multiple of 3 nucs in each to avoid frameshift 
# and checks that ALL species are represented
# then outputs a phylip file and fasta file with species and their full exon seqs
# python3 get_phylip_all_exons.py -i ./aligned_exons/ -o seqfile_pre-aligned
# NOTE for paml, you need to add some space between name and start of sequence (2 spaces), and add "I" for interlaved on top line
# NOTE for paml, names in trees will be wrong because phylip format cuts them off at 10 chars
# on server also

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
from Bio import SeqIO

import os
import sys
import argparse
import csv

def readFasta(fileName): 
    ''' helper
    '''
    recordsL = []
    for record in SeqIO.parse(fileName, "fasta"):
        recordsL.append(record)
    return recordsL

def readSpeciesL(speciesLFile):
    speciesL=[]
    with open(speciesLFile, 'r') as file:
        reader = csv.reader(file, delimiter=',')
        for row in reader:
            speciesL.append(row[0])
    return speciesL

def getRecordsTable(fastaNamesL, speciesL):
    ''' Takes a list of fasta names which each have a homologous exon, in order
    '''
    recordsTable = [speciesL]
    for fastaName in fastaNamesL:
        recordsL = readFasta(fastaName)
        if len(recordsL) != len(speciesL):
            # print(fastaName, " has only ", len(recordsL), " sequences out of ", len(speciesL), " species total")
            pass
        if len(recordsL) == len(speciesL): #if all species represented
            if validRecordsL(recordsL): # all sequences mult of 3 and all same len
                # add these sequences to their respective species lists
                print(fastaName)
                recordsTable.append(recordsL)
    return recordsTable


def makeFullRecords(recordsTable):
    ''' Takes records table and returns a list of records that have the full sequence for each species
    '''
    speciesL = recordsTable[0]
    sequenceL = ['' for i in range(len(speciesL))]
    for recordIndex in range(1,len(recordsTable)):
        for speciesIndex in range(len(speciesL)):
            thisSpecies = speciesL[speciesIndex]
            thisRecord = recordsTable[recordIndex][speciesIndex]
            thisSeq = str(thisRecord.seq)
            sequenceL[speciesIndex] = sequenceL[speciesIndex]+thisSeq
    # write recordsL
    fullRecordsL = []
    for i in range(len(speciesL)):
        thisSpecies = speciesL[i]
        thisSpeciesSeq = Seq(sequenceL[i], generic_dna)
        thisSpeciesRecord = SeqRecord(thisSpeciesSeq, id=thisSpecies, description='')
        fullRecordsL.append(thisSpeciesRecord)
    return fullRecordsL

def validRecordsL(recordsL):
    '''helper
    checks that for all records in recordsL, they are the same length which is a multiple of 3
    '''
    for record in recordsL:
        if len(record.seq)%3 != 0:
            print("not mult of 3")
            return False
    #print(record.name, " :this exon is valid")
    return True

def validRecordsLlen(recordsL):
    '''helper
    checks that for all records in recordsL, they are the same length which is a multiple of 3
    '''
    seqLens = []
    for record in recordsL:
        if len(record.seq)%3 != 0:
            print("not mult of 3")
            return False
        seqLens += [len(record.seq)]
    for i in range(len(seqLens)-1):
        # check all records are the same length
        if seqLens[i] != seqLens[i+1]:
            print("some species have difft length exon")
            return False
    print("this exon is valid")
    return True

def writePhylip(fullRecordsL, outname):
    # final fcn
    count = SeqIO.write(fullRecordsL, outname, "phylip")
    return

def writeFasta(fullRecordsL, outname):
    # final fcn
    count = SeqIO.write(fullRecordsL, outname, "fasta")
    return


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                   help="folder name",
                   type=str, default = "./homologous_fasta/")
    ap.add_argument("-s", "--speciesL",
            help="file name for output",
            type=str, default = "speciesL.txt")
    ap.add_argument("-o", "--outfile",
                help="file name for output",
                type=str, default = "seqfile")

    args = ap.parse_args()
    folderName = args.input
    outfile = args.outfile

    fastaNamesL = os.listdir(folderName)
    fastaNamesL = [folderName+fastaNamesL[i] for i in range(len(fastaNamesL))]
    # get the first exon first!
    fastaNamesL.reverse()
    speciesL = readSpeciesL(args.speciesL)
    recordsTable = getRecordsTable(fastaNamesL,speciesL)
    fullRecordsL = makeFullRecords(recordsTable)
    #writePhylip(fullRecordsL, outfile+".phylip")
    writeFasta(fullRecordsL, outfile+".fasta")
    

if __name__=="__main__":
    main()