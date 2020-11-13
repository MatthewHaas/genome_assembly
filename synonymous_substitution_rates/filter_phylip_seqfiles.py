import os
import sys
import re

infile = open(sys.argv[1], 'r') # seqfile in phylip format (but multiple lines per species)
outfile = open(sys.argv[2], 'w') # seqfile in phylip format (but one line per species)

seqFile = infile.readlines()

species = list() # initialize list to store species names, to keep track of those already seen

firstline = True
for line in seqFile:
        if firstline:
                outfile.write(line[0] + '  ' + line[2])
                firstline = False
                continue
        line.rstrip()
        line = line.split(' ')
        print(line[0])
        if line[0] in species:
                print(line[0])
                continue
        else:
                if line[0] == 'Osativa':
                        species.append(line[0])
                        print(line)
                        for item in line:
                                outfile.write(item)
                elif line[0] == 'Zlatifolia':
                        species.append(line[0])
                        for item in line:
                                outfile.write(item)
                elif line[0] == 'Zpalustris':
                        species.append(line[0])
                        for item in line:
                                outfile.write(item)
outfile.close()