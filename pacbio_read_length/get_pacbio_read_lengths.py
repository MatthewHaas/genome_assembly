import sys
import os
import matplotlib.pyplot as plt

infile = open(sys.argv[1], 'r')

data = infile.readlines()

read_lengths = list()

for line in data:
        line = line.split('_')
        length = int(line[3])
        read_lengths.append(length)


plt.hist(read_lengths, color = 'blue', edgecolor = 'black', bins = int(500))
plt.xlabel('PacBio Read Length (bp)')
plt.ylabel('Frequency')
plt.savefig('pacbio_length_distr.png')
