from Bio import SeqIO
import sys
import os

infile = open(sys.argv[1], 'r') # 'NWR_unique_orthogroup_list.txt'
outfile = open(sys.argv[2], 'w')

path = '/home/jkimball/shared/WR_Annotation/Orthofinder_rice_relatives/Grass_Proteomes/OrthoFinder/Results_Dec18/WorkingDirectory/OrthoFinder/Results_Dec21/Orthogroup_Sequences/'
orthogroup_list = infile.readlines()

for orthogroup in orthogroup_list:
	orthogroup = orthogroup.rstrip()
	for record in SeqIO.parse(str(path) + str(orthogroup)+'.fa', "fasta"):
	    outfile.write(record.id + '\n')
outfile.close()
