# Author: Fiona Callahan

# read in a fasta file that is input and output a list of coordinates where it finds gaps (ie N's)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein

import argparse
import os
import csv

def readSeqs(seqFileName):
    ''' Reads fasta file into string sequence and puts it in a dictionary with fasta ID as the key
    '''
    seqD = {}
    # Create nucleotide sequence objects
    for thisSeq in SeqIO.parse(seqFileName,"fasta"): # read the titin nucleotide sequence, will be stored as SeqRecord
        name = thisSeq.description
        mySeq = str(thisSeq.seq)
        seqD[name] = mySeq
    return seqD

"""
def readSeq(seqFileName):
    ''' Reads fasta file into a string sequence
    '''
    # Create nucleotide sequence objects
    ref_seq = SeqIO.read(seqFileName,"fasta") # read the titin nucleotide sequence, will be stored as SeqRecord
    mySeq = str(ref_seq.seq)
    return mySeq
"""

def find_gap_coords(mySeqsD):
    ''' input: dict of sequence -- values are sequences (strings)
        Output: dictionary of lists of coordinates of gap open and closes -- keyed by the same things as input
    '''
    gapsD = {}
    for name in mySeqsD.keys():
        mySeq = mySeqsD[name]
        gapL = []
        inGap = False
        for i in range(len(mySeq)):
            if mySeq[i] == "N" and inGap == False:
                # open gap
                inGap = True
                gapL.append([i,0])
            if mySeq[i] != "N" and inGap == True:
                # close gap
                inGap = False
                gapL[-1][1] = i
        gapsD[name] = gapL
    return gapsD

def countGapsStats(gapD, seqD):
    '''Takes gapD and seqD, two dictionaries with the same keys, and prints info about the percent of N's and the number of gaps
    '''
    nameL = []
    numGapsL = []
    numGapNucsL = []
    seqLenL = []
    for name in seqD.keys():
        print("Showing stats for: ", name)
        nameL.append(name)

        gapL = gapD[name]
        seqLen = len(seqD[name])

        gapLlen = len(gapL)
        print("The number of gaps is: ", gapLlen)
        numGapsL.append(gapLlen)

        numNucs = 0
        for gap in gapL:
            if(gap[1] == 0):
                gap[1] = seqLen
            numNucs += gap[1] - gap[0]
        print("Coordinates of the gaps are:", gapL)
        print("The number of gap nucleotides is: ", numNucs)
        numGapNucsL.append(numNucs)
        print("The total number of nucleotides is: ", seqLen)
        seqLenL.append(seqLen)
        percentGap = float(numNucs)/float(seqLen)
        print("For a proportion gap of: ", percentGap)
        print()

    # THIS CODE WORKS BUT THE IDEA DIDNT BECAUSE WHAT I DOWNLOADED ENDED UP BEING MOSTLY NOT GENOMIC
    
    outpath = "./"
    gapStats_outdir = outpath + "gapstats.csv"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    with open(gapStats_outdir, 'w') as csv_file:
        wr = csv.writer(csv_file, delimiter=',')
        wr.writerow(["name"]+nameL)
        wr.writerow(["numGaps"]+numGapsL)
        wr.writerow(["NumGapNucs"]+numGapNucsL)
        wr.writerow(["SeqLen"]+seqLenL)
        percentGapL = [numGapNucsL[i]/seqLenL[i] for i in range(len(seqLenL))]
        wr.writerow(["PercentGap"]+percentGapL)
    

def main():
    # get filename
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                   help="Path to the file(s) that contain the full titin nucleotide sequence",
                   type=str)

    args = ap.parse_args()
    fileName = args.input

    mySeqsD = readSeqs(fileName)
    gapsD = find_gap_coords(mySeqsD)
    countGapsStats(gapsD, mySeqsD)


if __name__=="__main__":
    main()