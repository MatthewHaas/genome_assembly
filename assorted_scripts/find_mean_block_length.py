import os
import sys
import numpy as np

infile = open(sys.argv[1], 'r')

data = infile.readlines()

lengths = list()
firstline = True
for line in data:
	if firstline:
		firstline = False
		continue
	line = line.strip()
	line = line.split('\t')
	block = float(line[2])
	lengths.append(block)

# Find mean of list
meanBlock = np.mean(lengths)

print('mean block length: ' + str(meanBlock))
