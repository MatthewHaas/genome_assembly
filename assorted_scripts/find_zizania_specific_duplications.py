import os
import sys

infile = open(sys.argv[1], 'r') # Duplications.tsv file (created by OrthoFinder)
outfile = open(sys.argv[2], 'w') # Give a name to write results to.

data = infile.readlines()

# Create empty dictionaries to hold data
genes1 = {}
genes2 = {}

firstline = True
for line in data:
    if firstline:
        firstline = False
        continue
    line.strip()
    line = line.split("\t")
    gene1 = line[5].split(',')
    gene2 = line[6].split(',')
    if line[1] == "N10":
        if len(gene1) == 2: 
            if len(gene2) == 2:
        #genes1[line[0]] = line[5]
        #genes2[line[0]] = line[6]
                outfile.write('\t'.join(line[0:])) # newline character not needed because it is already encoded in the input data (and was not removed using .strip()...)
        else:
            pass
outfile.close()
#print(genes1)
#print(genes2)
