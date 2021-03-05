from Bio import SeqIO
import os
import sys

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

for record in SeqIO.parse(infile, 'fasta'):
    if 'FUN' in record.description:
        outfile.write(record.format('fasta'))
        #print(record.format('fasta'))
    elif "Zlat" in record.description:
        outfile.write(record.format('fasta'))
        #print(record.format('fasta'))
    elif 'Os' in record.description:
        outfile.write(record.format('fasta'))
        #print(record.format('fasta'))

outfile.close()
