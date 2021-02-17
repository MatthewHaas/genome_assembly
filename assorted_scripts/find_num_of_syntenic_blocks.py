import os
import sys

infile = open(sys.argv[1], 'r') #open('201013_NWR_block_lengths.txt', 'r')

data = infile.readlines()

gene_counts = list()
firstline = True
for line in data:
	if firstline:
		firstline = False
		continue
	line = line.strip()
	line = line.split('\t')
	diff = line[2]
	gene_counts.append(float(diff))

total = 0
for item in gene_counts:
	total = total + item

print('total number of gene (pairs) is: ' + str(total))
