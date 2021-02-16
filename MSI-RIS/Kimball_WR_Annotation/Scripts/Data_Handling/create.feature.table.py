import sys

infile = open(sys.argv[1],'r') #PASA gff3 reordered
infile2 = open('blast2go_annot_WR.annot','r')#updated B2G annotations
infile3 = open('annotation_mapping_table_NCBI.txt','r')
infile4 = open('Zp_CDS_Classifications.txt','r')
infile5 = open('rice_gene_UTRs.txt','r')
infile6 = open('Zp_mRNA_Exon_Positions.txt','r')
outfile = open('feature.table.test.tbl','w')

mdict={}# map tid to gids
for line in infile3:
	f = line.strip().split('\t')
	ogid = f[2]
	otid = f[3]
	mdict[otid] = ogid

# saves blast2GO gene description information in dictionary. Key = gid, value = B2G description
b2g={}
for line in infile2:
	f = line.strip().split('\t')
	gid = f[0]
	if len(f) == 3:
		annot = f[2]
		b2g[gid] = annot

#print b2g

# save UTR info for genes
p3utrs = {}
p5utrs = {}
for line in infile5:
	f = line.strip().split('\t')
	gid = f[0]
	tid = f[1]
	utr5 = f[2]
	utr3 = f[3]
	if utr3 == '0':
		if gid in p3utrs:
			p3utrs[gid].append(tid)
		elif gid not in p3utrs:
			p3utrs[gid] = []
			p3utrs[gid].append(tid)
	if utr5 == '0':
		if gid in p5utrs:
			p5utrs[gid].append(tid)
		elif gid not in p5utrs:
			p5utrs[gid] = []
			p5utrs[gid].append(tid)

#print p5utrs

# save Tom's mRNA-exon data
mRNAs = {} # key = mRNA; values = list of exon start, stops
infile6.readline()
for line in infile6:
	f = line.strip().split('\t')
	tid = f[0]
	gid = f[1]
	exstarts = f[5].split(',')
	exstops = f[6].split(',')
	strand = f[7]
	mRNAs[tid] = []
	#print exstarts
	#print len(exstarts)
	for i in range(0,len(exstarts)): 
		#print i
		#print exstarts[i]
		#if gid/mRNA is on the positive strand and this is the first exon of the transcript
		if i == 0:
			#exstarts[i] = str(int(exstarts[i])+1)
			if gid in p5utrs:
				if tid in p5utrs[gid]:
					exstarts[i] = '<'+exstarts[i]
			if len(exstarts) == 1:
				if gid in p3utrs:
					if tid in p3utrs[gid]:
						exstops[i] = '>'+exstops[i]
		elif i == (len(exstarts)-1):
			if gid in p3utrs:
				if tid in p3utrs[gid]:
					exstops[i] = '>'+exstops[i]	
		mRNAs[tid].append([exstarts[i],exstops[i]])
#print mRNAs



# load Tom's CDS data#infile4.readline()
infile4.readline()
CDSdict={}
cds_phase = {}
for line in infile4:
	f =line.strip().split('\t')
	cdsid = f[0]# CDS ID
	starts = f[6].split(',')
	stops = f[7].split(',')
	CDSdict[cdsid] = []
	for i in range(0,len(starts)):
		if i == 0 and cdsid not in cds_phase:
			phasef = f[9].split(',')
			phase = phasef[0]
			if phase != '0':
				phase = str(int(phase) +1)# removed + 1
				cds_phase[cdsid] = phase
		"""# check if symbol in start/stop
		t_up = 0
		t_down = 0
		if '<' in starts[i]:
			print starts[i]
			starts[i] =starts[i][1:]
			print starts[i]
			t_up = 1
		if '>' in stops[i]:
			print stops[i]
			stops[i] = stops[i][1:]
			print stops[i]
			t_down = 1
		print starts[i],stops[i]
		if int(starts[i]) > int(stops[i]):
			stops[i] = str(int(stops[i]) + 1)
			if t_down == 1:
				stops[i] = '>'+ stops[i]
			if t_up == 1:
				starts[i] == '<'+str(starts[i])
		elif int(starts[i]) < int(stops[i]):
			starts[i] = str(int(starts[i]) + 1)
			if t_down == 1:
				stops[i] = '>'+ str(stops[i])
			if t_up == 1:
				starts[i] = '<'+ starts[i]
		"""
		CDSdict[cdsid].append([starts[i],stops[i]])

#print CDSdict
print cds_phase
print len(cds_phase)

for cds in cds_phase:
	print cds

CDStrack = {}
chromdict = {}
for line in infile:
	#print line
	if line[0] != '#' and line[0] != '' and line[0] != '\n':
		f = line.strip().split('\t')
		#print f
		chrom = f[0]
		feat = f[2]
		start = f[3]
		stop  = f[4]
		strand = f[6]
		phase = f[7]
		if chrom not in chromdict:
			chromdict[chrom] = []
			outfile.write('>Feature\t'+chrom+'\n')
			if strand == '-': # reverse the start and stop of the feature
				start = f[4]
				stop = f[3]
				if feat == 'gene':
					f2 = f[8].split(';')
					gid = f2[0][3:]
					if gid in p3utrs:
						stop = '>'+stop
					if gid in p5utrs:
						start = '<'+start	 
					outfile.write(start+'\t'+stop+'\t'+feat+'\n')
					outfile.write('\t\t\t'+'locus_tag'+'\t'+gid+'\n')
					
			elif strand != '-': # start and stop of features are not reversed
				if feat == 'gene':
					f2 = f[8].split(';')
					gid = f2[0][3:]
					if gid in p3utrs:
						stop = '>'+stop 
					if gid in p5utrs:
						start = '<'+start
					outfile.write(start+'\t'+stop+'\t'+feat+'\n')
					outfile.write('\t\t\t'+'locus_tag'+'\t'+gid+'\n')
					
		elif chrom in chromdict:
			if strand == '-': # reverse the start and stop of the feature
				start = f[4]
				stop = f[3]
				if feat == 'gene':
					f2 = f[8].split(';')
					gid = f2[0][3:]
					if gid in p3utrs:
                                                stop = '>'+stop
                                        if gid in p5utrs:
                                                start = '<'+start
					outfile.write(start+'\t'+stop+'\t'+feat+'\n')
					outfile.write('\t\t\t'+'locus_tag'+'\t'+gid+'\n')
				elif feat == 'mRNA':
					f2 = f[8].split(';')
					tid = f2[0][3:]
					gid = mdict[tid]
					# remove the follow five lines to make the first exon of the mRNA the first entry, instead of the mRNA's start and stop
					#if gid in p3utrs and tid in p3utrs[gid]:
					#	stop = '>'+stop
					#if gid in p5utrs and tid in p5utrs[gid]:
					#	start = '<'+start	
					#outfile.write(start+'\t'+stop+'\t'+feat+'\n')
					exons = mRNAs[tid]
					etmp = 0
					for exon in exons:
						etmp += 1
						[estart,estop] = exon
						if etmp == 1:
							outfile.write(estart+'\t'+estop+'\tmRNA\n')
						else:
							outfile.write(estart+'\t'+estop+'\n')
					outfile.write('\t\t\t'+'transcript_id'+'\t'+tid+'\n')
					outfile.write('\t\t\t'+'protein_id'+'\t'+tid+'\n')
					if gid in b2g:
						outfile.write('\t\t\t'+'product'+'\t'+b2g[gid]+'\n')
					elif gid not in b2g:
						outfile.write('\t\t\t'+'product'+'\t'+'predicted protein'+'\n')
				elif feat == 'CDS':
					f2 = f[8].split(';')
					tid = f2[1][7:]	
					gid = mdict[tid]
					cdsID = f2[0][3:] 
					#print cdsID
					if cdsID not in CDStrack:
						CDStrack[cdsID] = 0
 						CDSes = CDSdict[cdsID]
						#print CDSes
						for i in range(0,len(CDSes)):
							[cdsstart,cdsstop] = CDSes[i]
							if i == 0:
								outfile.write(cdsstart+'\t'+cdsstop+'\tCDS\n')
							else:
								outfile.write(cdsstart+'\t'+cdsstop+'\n')
						#if '3_prime_partial' in line:
						#	outfile.write(start+'\t>'+stop+'\t'+feat+'\n')
						#elif '5_prime_partial' in line:
						#	outfile.write('<'+start+'\t'+stop+'\t'+feat+'\n')
						#elif '3_prime_partial' not in line and '5_prime_partial' not in line:
						#	outfile.write(start+'\t'+stop+'\t'+feat+'\n')
						if cdsID in cds_phase:
							outfile.write('\t\t\tcodon_start\t'+cds_phase[cdsID]+'\n')
						outfile.write('\t\t\t'+'transcript_id'+'\t'+tid+'\n')
						outfile.write('\t\t\t'+'protein_id'+'\t'+tid+'\n')
						#if int(phase) in [1,2]:
						#	outfile.write('\t\t\t'+'codon_start'+'\t'+str(phase)+'\n')
						if gid in b2g:
							outfile.write('\t\t\t'+'product'+'\t'+b2g[gid]+'\n')
						elif gid not in b2g:
							outfile.write('\t\t\t'+'product'+'\t'+'predicted protein'+'\n')
			elif strand != '-': # start and stop of features are not reversed
				if feat == 'gene':
					f2 = f[8].split(';')
					gid = f2[0][3:]
					if gid in p3utrs:
                                                stop = '>'+stop
                                        if gid in p5utrs:
                                                start = '<'+start
					outfile.write(start+'\t'+stop+'\t'+feat+'\n')
					outfile.write('\t\t\t'+'locus_tag'+'\t'+gid+'\n')
				elif feat == 'mRNA':
					f2 = f[8].split(';')
					tid = f2[0][3:]
					gid = mdict[tid]
					# remove the follow five lines to make the first exon of the mRNA the first entry, instead of the mRNA's start and stop
					#if gid in p3utrs and tid in p3utrs[gid]:
					#	stop = '>'+stop
					#if gid in p5utrs and tid in p5utrs[gid]:
					#	start = '<'+start
					#outfile.write(start+'\t'+stop+'\t'+feat+'\n')
					exons = mRNAs[tid]
					etmp=0
					for exon in exons:
						etmp += 1
						[estart,estop] = exon
						if etmp == 1:
							outfile.write(estart+'\t'+estop+'\tmRNA\n')
						else:
							outfile.write(estart+'\t'+estop+'\n')
					outfile.write('\t\t\t'+'transcript_id'+'\t'+tid+'\n')
					outfile.write('\t\t\t'+'protein_id'+'\t'+tid+'\n')
					if gid in b2g:
						outfile.write('\t\t\t'+'product'+'\t'+b2g[gid]+'\n')
					elif gid not in b2g:
						outfile.write('\t\t\t'+'product'+'\t'+'predicted protein'+'\n')
				elif feat == 'CDS':
					f2 = f[8].split(';')
					f3 = f2[0].split('.')
					tid = f2[1][7:]
					gid = mdict[tid]
					cdsID = f2[0][3:]
					#print cdsID
					if cdsID not in CDStrack:
						CDStrack[cdsID] = 0
						CDSes = CDSdict[cdsID]
						#print CDSes
						for i in range(0,len(CDSes)):
							[cdsstart,cdsstop] = CDSes[i]
							if i == 0:
								outfile.write(cdsstart+'\t'+cdsstop+'\tCDS\n')
							else:
								outfile.write(cdsstart+'\t'+cdsstop+'\n')

					#if '3_prime_partial' in line:
					#	outfile.write(start+'\t>'+stop+'\t'+feat+'\n')
					#elif '5_prime_partial' in line:
					#	outfile.write('<'+start+'\t'+stop+'\t'+feat+'\n')
					#elif '3_prime_partial' not in line and '5_prime_partial' not in line:
					#	outfile.write(start+'\t'+stop+'\t'+feat+'\n')
						if cdsID in cds_phase:
							outfile.write('\t\t\tcodon_start\t'+cds_phase[cdsID]+'\n')	
						outfile.write('\t\t\t'+'transcript_id'+'\t'+tid+'\n')
						outfile.write('\t\t\t'+'protein_id'+'\t'+tid+'\n')
						#if int(phase) in [1,2]:
						#	outfile.write('\t\t\t'+'codon_start'+'\t'+str(phase)+'\n')
						if gid in b2g:
							outfile.write('\t\t\t'+'product'+'\t'+b2g[gid]+'\n')
						elif gid not in b2g:
							outfile.write('\t\t\t'+'product'+'\t'+'predicted protein'+'\n')
