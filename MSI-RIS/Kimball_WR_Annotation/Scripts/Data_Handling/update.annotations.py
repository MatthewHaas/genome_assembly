import sys
# updates the feature table and gff3 file for NCBI

infile = open('rice.gene_structures_post_PASA_updates.21917.gff3','r')
infile2= open('../genome/contig_mapping_table2.txt','r')
infile3 = open('feature.table.test.tbl','r')
outfile = open('pasa_gene_models_NCBI.gff3','w')
outfile2 = open('feature.table.NCBI.tbl','w')
outfile3 = open('annotation_mapping_table_NCBI.txt','w')

scafs = {}
for line in infile2:
	f = line.strip().split('\t')
	oid = f[0]
	nid = f[1]
	scafs[oid] = nid


gdict = {} # keys are old gIDs, values are new geneIDs
gtdict = {} # keys are old gIDs, values old transcript IDs
counter = 0
for line in infile:
	if line[0] != '#' and line[0] != '\n':
		f = line.strip().split('\t')
		chrom = f[0]
		nchrom = scafs[chrom]
		feat = f[2]
		if feat == 'gene':
			counter += 1
			f2 = f[8].split(';')
			gid = f2[0][3:]
			ngid = nchrom+'g'+str(counter)
			gdict[gid] = ngid
		elif feat == 'mRNA':
			f2 = f[8].split(';')
			tid = f2[0][3:]
			gid = f2[1][7:]
			if gid in gtdict:
				gtdict[gid].append(tid)
			elif gid not in gtdict:
				gtdict[gid] = []
				gtdict[gid].append(tid)
			

tdict={} # old transcript IDs [keys] and new-transcript IDs[values]
for gid in gtdict:
	trans = gtdict[gid]
	ngid = gdict[gid]
	c = 0
	for t in trans:
		c += 1
		tdict[t] = ngid+'-'+str(c)
		outfile3.write(ngid+'\t'+ngid+'-'+str(c)+'\t'+gid+'\t'+t+'\n')
		
	
for line in infile3:
	if line[:8] == '>Feature':
		f = line.strip().split('\t')
		ochr = f[1]
		nchr = scafs[ochr]
		outfile2.write(f[0]+'\t'+nchr+'\n')
	else:
		f = line.strip().split('\t')
		#print f
		if 'transcript_id' in line or 'protein_id' in line or 'locus_tag' in line:	
		#if len(f) == 2:
			if f[0] == "transcript_id":
				outfile2.write('\t\t\ttranscript_id\tgnl|KimballUMN|'+tdict[f[1]]+'\n')
			elif f[0] == "protein_id":
				outfile2.write('\t\t\tprotein_id\tgnl|KimballUMN|'+tdict[f[1]]+'\n')
			elif f[0] == "locus_tag":
				outfile2.write('\t\t\tlocus_tag\t'+gdict[f[1]]+'\n')
			else:
				outfile2.write(line)
		else:
			outfile2.write(line)
