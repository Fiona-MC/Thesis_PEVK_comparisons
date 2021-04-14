# Author: Fiona Callahan
# get pevk stats
# want to know how many exons and how long they are total
# use output of PEVK_finder CSVs
# GOES THROUGH WHOLE FOLDER TO COMPILE THE DATA

# Note: if you run it the second time it wont work because it will try to run the csv file that it deposits in the folder with the input


import csv
import argparse
import os

def readfile(filename):
    # read CSV file in:
    # rows are names, lengths, ratios
    print(filename)
    with open(filename, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        data = list(reader)
    return data

def getDataList(data, speciesName):
    '''input: species name and the output of reading the CSV with exon lengths and ratios in
    output: list with species name, number of exons, and total exon length
    '''
    exonLens = data[1]
    numExons = 0
    totalexonLen = 0
    for exon in exonLens:
        numExons += 1
        totalexonLen += int(exon)
    return [speciesName, numExons, totalexonLen]

def writeCsvRow(listData,outfile):
    """writes new line in CSV file with the list of things in it: in this case [speciesName, numExons, totalexonLen]
    """
    with open(outfile, 'a') as csv_file:
        wr = csv.writer(csv_file, delimiter=',')
        wr.writerow(listData)

def main():
    # get filename
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                   help="Path to the folder that contains the exon lengths and ratios CSVs you want -- output of pevk_finder with extra outputs",
                   type=str)
    ap.add_argument("-l", "--label",
                   help="label for file",
                   type=str)
    ap.add_argument("-o", "--outdir",
                   help="path to the place you want the file",
                   type=str)             

    args = ap.parse_args()
    for filename in os.listdir(args.input):
        if filename.endswith('.csv'):
            filepath = args.input + filename

            speciesName = filename # really I should change this
            outfile = args.outdir+args.label+"_pevk_stats.csv"

            data = readfile(filepath)
            dataList = getDataList(data, speciesName)
            writeCsvRow(dataList,outfile)

if __name__=="__main__":
    main()
