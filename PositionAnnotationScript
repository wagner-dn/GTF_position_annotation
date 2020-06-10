#!/usr/bin/env python3

# This program annotates a list of genome positions from an input file of "Chr<TAB>Position" format
# and a GTF file as reference to an output file in the current directory.
# Usage: ./position_annotation_gtf_DNW_faster.py CoordinateFile GTF_file
### Note, the output file will be placed in the same directory as the CoordinateFile
# Run by moving script to desired directory and running
# or calling absolute path to script from desired directory.

### Written by Dominique Neitzel Wagner, 29. October 2019
### modified by DNW, 11. November 2019

import sys
import pandas as pd
# from datetime import datetime

#Initialize timing
# startTime = datetime.now()

COORDFileName = sys.argv[1]
GTFFileName = sys.argv[2]

# COORDFileName = 'sample_files/annotate/coordinates_to_annotate.txt'
# GTFFileName = 'sample_files/gtf/hg19_annotations.gtf'

OutFileName = COORDFileName[:-4] + "-annotated_faster_v5.txt"
OutFile=open(OutFileName, 'w')

# Load data frame with panda because it's fast and can index easily
COORDFile = pd.read_csv(COORDFileName, header=None, sep='\t',
	names=['CHR', 'POS'])
COORDFile_sort = COORDFile.sort_values(['CHR','POS'])
COORDFile_sort = COORDFile_sort.reset_index(drop=True)

# Create dictionary of chromosomes to increase search speed
# First version didn't have this and was MUCH slower
GTFdictionary = {}
with open(GTFFileName) as infile:
	for line in infile:
		newline = line.strip().split("\t")
		chrom = newline[0]
		if chrom not in GTFdictionary:
			GTFdictionary[chrom] = []
		newline[3] = int(newline[3])
		newline[4] = int(newline[4])
		GTFdictionary[chrom].append(newline)

# Initialize output dictionary
outputDictionary = {}
COORDLineNumber = 0
GTFLineNumber = 0

# Loop through Coordinate file for chromosome dictionary IDs and find matches in range from GTF file
for COORDLineNumber in range(len(COORDFile_sort)):
	RefCHR = COORDFile_sort.loc[COORDLineNumber, 'CHR']
	RefPOS = COORDFile_sort.loc[COORDLineNumber, 'POS']
	if (RefCHR,RefPOS) not in outputDictionary: # Initialize new positions
		outputDictionary[(RefCHR,RefPOS)] = "No Match Found" # Default output
	if RefCHR in GTFdictionary:
		GTFlines_currentChrom = GTFdictionary[RefCHR]
		for GTFLineNumber in range(len(GTFlines_currentChrom)):
			if (GTFlines_currentChrom[GTFLineNumber][3] <= RefPOS) and (GTFlines_currentChrom[GTFLineNumber][4]+1 >= RefPOS):
				GeneAttributes = GTFlines_currentChrom[GTFLineNumber][8]
				geneID = GeneAttributes.split(";")[0].split()[1][1:-1]
				outputDictionary[(RefCHR,RefPOS)] = geneID # Override default with geneID
for pos in outputDictionary: # Iterate over all of the keys in the entire dictionary
	output = "\t".join(map(str,pos)) # Convert the tuple used as a key into a joined string
	OutFile.write(output+"\t"+str(outputDictionary[pos])+"\n") # tack on the dictionary value, and line break

# Close output file (or else it will be empty!)
OutFile.close()
# print(datetime.now() - startTime)
