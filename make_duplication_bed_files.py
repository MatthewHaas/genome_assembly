import os

# Lots of files to open here:
# First, we need to open the original BED files.
# Lastly, we need to open new BED files to write output to.

infile1 = open('wild_rice.bed', 'r')
infile2 = open('latifolia.bed', 'r')

infile3 = open('duplications_gene1.csv', 'r')
infile4 = open('duplications_gene2.csv', 'r')

# Output files. We want one per species and gene1/gene2 combination
wildRiceDuplicationsGene1Bed = open('wild_rice_duplications_gene1.bed', 'w') 
wildRiceDuplicationsGene2Bed = open('wild_rice_duplications_gene2.bed', 'w')
latifoliaDuplicationsGene1Bed = open('latifolia_duplications_gene1.bed', 'w')
latifoliaDuplicationsGene2Bed = open('latifolia_duplications_gene2.bed', 'w')

wildRiceBed = infile1.readlines()
latifoliaBed = infile2.readlines()

gene1_duplications = infile3.readlines()
gene2_duplications = infile4.readlines()

# Initialize lists
gene1_species1 = []
gene1_species2 = []
gene2_species1 = []
gene2_species2 = []

for line in gene1_duplications:
	line.rstrip()
	line = line.split(',')
	print(line)
	gene1_species1.append(line[0])
	gene1_species2.append(line[1])
	
for line in gene2_duplications:
	line.rstrip()
	line = line.split(',')
	gene2_species1.append(line[0])
	gene2_species2.append(line[1])
	


for line in wildRiceBed:
	line.rstrip()
	line = line.split('\t')
	if line[3] in gene1_species1:
		wildRiceDuplicationsGene1Bed.write('\t'.join(line[0:]))
	elif line[3] in gene1_species2:
		wildRiceDuplicationsGene1Bed.write('\t'.join(line[0:]))
	elif line[3] in gene2_species1:
		wildRiceDuplicationsGene2Bed.write('\t'.join(line[0:]))
	elif line[3] in gene2_species2:
		wildRiceDuplicationsGene2Bed.write('\t'.join(line[0:]))
		
for line in latifoliaBed:
	line.rstrip()
	line = line.split('\t')
	if line[3] in gene1_species1:
		latifoliaDuplicationsGene1Bed.write('\t'.join(line[0:]))
	elif line[3] in gene1_species2:
		latifoliaDuplicationsGene2Bed.write('\t'.join(line[0:]))
