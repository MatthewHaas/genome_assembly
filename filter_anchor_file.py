import os
import sys

infile = open(sys.argv[1], 'r') # this should be the anchor file
infile2 = open(sys.argv[2], 'r') # this should be the file containing gene/chr info, to make a dictionary (/Users/matthewwilliamhaas/Documents/wild_rice/syntenic_genes.csv)
outfile = open(sys.argv[3], 'w') # this should be the new anchor file

anchors = infile.readlines()
genes = infile2.readlines()

# Makes a list of all genes from annotation file that are on chr1-15 plus scaffold 16 and scaffold 458
gene_list = list()
for line in genes:
	gene = line.rstrip()
	gene_list.append(gene)

outfile.write('###' + '\n') # put the header back in
firstline = True
for line in anchors:
	if firstline:
		firstline = False
		continue
	line = line.strip()
	line = line.split('\t')
	if line[0] in gene_list:
		outfile.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\n')
	else:
		continue
outfile.close()
