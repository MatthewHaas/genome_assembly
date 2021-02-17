import os
import sys

infile = open('./tab-separated_value_files/Duplications_full.tsv', 'r')
infile2 = open('starting_gene_in_blocks.txt', 'r')
outfile = open('duplicated_blocks.tsv', 'w')

Duplications = infile.readlines()
StartOfGeneBlocks = infile2.readlines() #['Zizania_palustris_FUN_018861-T1', 'Zizania_palustris_FUN_001795-T1']

outfile.write('Orthogroup' + '\t' + 'Genes 1' + '\t' + 'Genes 2' + '\n')

for gene in StartOfGeneBlocks:
    print(gene)
    firstline = True
    for line in Duplications:
        if firstline:
            firstline = False
            continue
        line.strip()
        line = line.split('\t')
        genes1 = line[5].split(',')
        genes2 = line[6].split(',')
        if gene in genes1:
            outfile.write(line[0] + '\t' + line[5] + '\t' + line[6] + '\n')
        elif gene in genes2:
            outfile.write(line[0] + '\t' + line[5] + '\t' + line[6] + '\n')

outfile.close()
