from Bio import SeqIO
import os
import sys


infile = open(sys.argv[1], 'r')
infile2 = open(sys.argv[2], 'r')
outfile = open(sys.argv[3], 'w')

duplication_file = infile2.readlines()

pairs = list()
for line in duplication_file:
        line = line.rstrip('\n')
        line = line.split(',')
        orthogroup = line[0]
        pairs.append((line[1], line[2]))
print(pairs)

for record in SeqIO.parse(infile, 'fasta'):
        for a,b in pairs:
                if a in record.description:
                        outfile.write(record.format('fasta'))
                elif b in record.description:
                        outfile.write(record.format('fasta'))

outfile.close()
