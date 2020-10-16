import os
import sys

infile = open(sys.argv[1], 'r') # *simple file from MCscan, no header
infile2 = open(sys.argv[2], 'r') # BED file, no header
outfile = open(sys.argv[3], 'w') # The file to which we will write the output (chromosome and block length)

blocks = infile.readlines()
bed_file = infile2.readlines()

outfile.write('chr' + '\t' + 'block_length' + '\n') # header
for gene_pair in blocks:
	gene_pair = gene_pair.split('\t')
	start_gene = gene_pair[0] # gene_pair[0] accesses first genome, gene_pair[2] accesses second genome
	end_gene = gene_pair[1] # gene_pair[1] accesses first genome, gene_pair[3] accesses second genome
	for line in bed_file:
		line = line.split('\t')
		if line[3] == start_gene:
			start_pos = float(line[1])
		elif line[3] == end_gene:
			end_pos = float(line[2])
		elif start_gene in line[3]:
			start_pos = float(line[1])
		elif end_gene in line[3]:
			end_pos = float(line[2])
		else:
			continue
		try:
			start_pos and end_pos
		except NameError:
			continue
		if end_pos > start_pos:
			block_length = end_pos - start_pos
			block_div_1000 = block_length/1000
			outfile.write(line[0] + '\t' + str(block_length) + '\t' + str(block_div_1000) + '\n')
		else:
			continue
