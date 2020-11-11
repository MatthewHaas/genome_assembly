import os
import sys
import re

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

file = infile.readlines()

firstline = True
for line in file:
	if firstline:
		outfile.write(line)
		firstline = False
		continue
	if 'Os' in line:
		line = re.sub('^Os.+[0-9]+', 'Osativa  ', line)
		outfile.write(line)
	elif 'Zlat' in line:
		line = re.sub('^Zlat_[0-9]+', 'Zlatifolia  ', line)
		outfile.write(line)
	elif 'FUN' in line:
		line = re.sub('^FUN_[0-9]+', 'Zpalustris  ', line)
		outfile.write(line)

outfile.close()
