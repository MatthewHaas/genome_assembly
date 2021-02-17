import os
import sys

infile = open(sys.argv[1], 'r') # input bed file

outfile = open(sys.argv[2], 'w') # output bed file

bedFile = infile.readlines()

majorChroms = ['ZPchr0001', 'ZPchr0002', 'ZPchr0003', 'ZPchr0004', 'ZPchr0005', 'ZPchr0006', 'ZPchr0007', 'ZPchr0008', 'ZPchr0009',
			  'ZPchr0010', 'ZPchr0011', 'ZPchr0012', 'ZPchr0013', 'ZPchr0014', 'ZPchr0015', 'ZPchr0016', 'ZPchr0458']

firstline = True
for line in bedFile:
	if firstline:
		firstline = False
		continue
	line = line.strip()
	line = line.split('\t')
	if line[0] in majorChroms:
		outfile.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\n')
	else:
		continue
outfile.close()
