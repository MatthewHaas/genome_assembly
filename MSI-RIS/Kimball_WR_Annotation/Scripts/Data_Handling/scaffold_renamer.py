import sys

infile = open('zizania_palustris_13Nov2018_okGsv_renamedNCBI.fasta','r')
infile2 = open('contig_mapping_table.txt','r')
outfile = open('zizania_palustris_13Nov2018_okGsv_renamedNCBI2.fasta','w')
outfile2 = open('contig_mapping_table2.txt','w')

for line in infile:
	if line[0] == '>':
		if len(line.strip()) == 5:
			cid = line[4:].strip()
			nid = '>ZPchr000'+cid+'\n'
		elif len(line.strip()) == 6:
			cid = line[4:].strip()
			nid = '>ZPchr00'+cid+'\n'
		elif len(line.strip()) == 7:
			cid = line[4:].strip()
			nid = '>ZPchr0'+cid+'\n'
		else:
			
			cid = line[4:].strip()
			nid = '>ZPchr'+cid+'\n'
		outfile.write(nid)
	else:
		outfile.write(line)

for line in infile2:
	f = line.strip().split('\t')
	cid = f[1][3:]
	if len(cid) == 1:
		nid = 'ZPchr000'+cid	
	elif len(cid) == 2:
		nid = 'ZPchr00'+cid
	elif len(cid) == 3:
		nid = 'ZPchr0'+cid
	elif len(cid) == 4:
		nid = 'ZPchr'+cid
	outfile2.write(f[0]+'\t'+nid+'\n')
