import re
import os

os.chdir('/home/jkimball/haasx092/other_synteny_figures/updated_dot_plot')

infile1 = 'rice.gene_structures_post_PASA_updates.21917.gff3'
infile2 = 'rice.gene_structures_post_PASA_updates.nucleotides.fsa'
gene_conv = 'contig_mapping_table2.txt'
outfile1 = 'rice.gene_structures_post_PASA_updates_ReformattedNames.gff3'

gene_name_dict = {}
with open(gene_conv) as gene_names:
	for line in gene_names:
		key, val = line.split()
		gene_name_dict[key] = val

input =  open(infile1, 'r')
output = open(outfile1, 'w')
data = input.readlines()
for line in data:
	line = line.split('\t')
	if line[0] in gene_name_dict:
		line[0] = line[0].replace(str(line[0]), str(gene_name_dict[line[0]]))
		listToStr = ' '.join(map(str, line))
		output.write(listToStr)
output.close()
