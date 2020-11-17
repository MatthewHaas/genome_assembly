from Bio import SeqIO
import os
import sys


infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

pairs = [('Zlat_10041350','FUN_011478-T1')]

for record in SeqIO.parse(infile, 'fasta'):
        for a,b in pairs:
                if a in record.description:
                        outfile.write(record.format('fasta'))
                elif b in record.description:
                        outfile.write(record.format('fasta'))

outfile.close()
