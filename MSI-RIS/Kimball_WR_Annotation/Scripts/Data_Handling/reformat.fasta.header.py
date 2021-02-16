import sys

infile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')

for line in infile:
	if line[0] == '>':
		f = line.strip().split(';')
		outfile.write(f[0]+'\n')
	elif line[0] != '>':
		outfile.write(line)
