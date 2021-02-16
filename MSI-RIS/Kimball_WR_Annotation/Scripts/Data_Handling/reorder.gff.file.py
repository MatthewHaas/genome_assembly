import sys

# reorder rice PASA GFF file to fix genes that are out of order on the scaffolds

infile=open('genes.sorted.bed','r')
infile2=open('rice.gene_structures_post_PASA_updates.21917.gff3','r')
outfile = open('rice.gene_structures_post_PASA_updates.21917.reordered.gff3','w')

gdict={}
seq = ''
for line in infile2:
	if line[0] != '#' and line[0] != '\n':
		f= line.strip().split('\t')
		feat = f[2]
		if feat == 'gene' and seq != '':
			gdict[gid] = seq
			gidf = f[8].split(';')
			gid = gidf[0]
			seq = line
		elif feat == 'gene' and seq == '':
			gidf = f[8].split(';')
			gid = gidf[0]
			seq += line	
		else:
			seq += line

gdict[gid] = seq

for line in infile:
	f = line.strip().split('\t')
	gid = f[3]
	outfile.write(gdict[gid])
