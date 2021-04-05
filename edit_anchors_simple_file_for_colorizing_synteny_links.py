import os
import sys

infile = open(sys.argv[1], 'r')

infile = open('wild_rice.oryza.anchors.simple', 'r')
infile2 = open('NWR_chr6_genes_old_names.txt', 'r')
outfile = open('wild_rice.oryza.anchors2.simple', 'w')


data = infile.readlines()
genes = infile2.readlines()

geneList = list()
for line in genes:
	line = line.strip()
	geneList.append(line)

for line in data:
	line = line.strip()
	line = line.split('\t')
	if line[0] in geneList:
		modified_gene1_field = 'r*' + str(line[0])
		outfile.write(modified_gene1_field + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\n')
	else:
		outfile.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\n')
		
outfile.close()
