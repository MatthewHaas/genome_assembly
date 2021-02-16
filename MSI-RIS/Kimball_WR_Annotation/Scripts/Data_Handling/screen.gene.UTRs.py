import sys

infile = open('rice.gene_structures_post_PASA_updates.21917.gff3','r')
infile2 = open('annotation_mapping_table_NCBI.txt','r')
outfile = open('rice_gene_UTRs.txt','w')

mdict = {}# keys are old transcripts; values are old genes
for line in infile2:
	f = line.strip().split('\t')
	gid = f[2]
	tid = f[3]
	mdict[tid] = gid

gdict = {}
for line in infile:
	if line[0] != '#' and line[0] != '\n':
		f =line.strip().split('\t')
		feat = f[2]
		if feat == 'gene':
			f2 = f[8].split(';')
			gid = f2[0][3:]
			gdict[gid] = {}
			#gdict[gid] = [0,0]
		elif feat == 'mRNA':
			f2 = f[8].split(';')
			tid = f2[0][3:]
			if gid in gdict and tid not in gdict[gid]:
				gdict[gid][tid] = [0,0]
		elif feat == 'five_prime_UTR':
			f2 = f[8].split(';')
			tid = f2[1][7:]
			gid = mdict[tid]	
			if gid in gdict and tid in gdict[gid]:
				gdict[gid][tid][0] += 1
		elif feat == 'three_prime_UTR':
			f2 = f[8].split(';')
			tid = f2[1][7:]
			gid = mdict[tid]
			if gid in gdict and tid in gdict[gid]:
				gdict[gid][tid][1]+=1

outfile.write('gid\ttid\t5_UTR\t3_UTR\n')
for g in gdict:
	tids = gdict[g]
	for t in tids:
		outfile.write(g+'\t'+t+'\t'+str(gdict[g][t][0])+'\t'+str(gdict[g][t][1])+'\n')
