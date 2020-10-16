import os
import re
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

infile1 = open(str(sys.argv[1]), 'r') #open('wild_rice.oryza.anchors.simple', 'r')
infile2 = open(str(sys.argv[2]), 'r') #open('wild_rice.bed', 'r')
outfile = open(str(sys.argv[3]), 'w') #open('NWR_blocks_start_end.txt', 'w')
outfile2 = open(str(sys.argv[4]), 'w')#open('201013_NWR_block_lengths.txt', 'w')

data = infile1.readlines()

outfile.write('start_gene' +  '\t' + 'end_gene' + '\n')
for line in data:
	line = line.split('\t')
	start_gene = line[2] # should be line[0] to access first genome or line[2] to access second genome
	end_gene = line[3] # should be line[1] to access first genome or line[3] to access second genome
	outfile.write(start_gene + '\t' +  end_gene + '\n')

outfile.close()

NWR_bed = infile2.readlines()	
outfile = open(sys.argv[3],'r') #open('NWR_blocks_start_end.txt', 'r')	
blocks = outfile.readlines()

outfile2.write('start'+ '\t' + 'end' + '\t' + 'diff' + '\n') 
firstline = True
for item in blocks:
	item = re.split('_|\t|\n', item)
	if firstline:
		firstline = False
		continue
	start = item[1]
	end = item[3]
	if 'g' in start:
		start = start[5:]
	if 'g' in end:
		end = end[5:]
	diff = float(end) - float(start)
	outfile2.write(start + '\t' +  end + '\t' + str(diff) + '\n')

outfile2.close()
# Make a plot showing the distribution of block sizes

infile3 = open(sys.argv[4], 'r') #open('201013_NWR_block_lengths.txt', 'r')
blockData = infile3.readlines()
blockSize = list()
firstline = True
for item in blockData:
	item = re.split('\t', item)
	if firstline:
		firstline = False
		continue
	val = float(item[2])
	blockSize.append(val)

mean = round(np.mean(blockSize), 2)
stdev = round(np.std(blockSize), 2)
stats = 'mean: ' + str(mean) + '\n' + 'stdev: ' + str(stdev)
	
#plt.plot()
#plt.hist(sorted(blockSize), bins = 50, color = '#235e39')
#plt.xlabel('Number of genes in block')
#plt.ylabel('Number of blocks')
#plt.text(x = 100, y = 100, s = stats)
#plt.savefig('2010013_genes_per_block_distribution_MCscan.png', dpi=300)
