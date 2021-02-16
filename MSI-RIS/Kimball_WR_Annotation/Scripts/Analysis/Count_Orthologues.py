#!/usr/bin/env python
"""Simple script to count the number of genes in each orthologous group that
are recorded in an orthologues.txt file from Orthofinder. Takes one argument:
    1) orthologues.txt (gzipped)"""

import sys
import gzip

try:
    orthogz = sys.argv[1]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)

with gzip.open(orthogz, 'rt') as f:
    for index, line in enumerate(f):
        if index == 0:
            header = line.strip().split('\t')
            print(','.join(header))
        else:
            tmp = line.strip().split('\t')
            ogname = tmp[0]
            og1_genes = tmp[1].split(', ')
            og2_genes = tmp[2].split(', ')
            og1_count = str(len(og1_genes))
            og2_count = str(len(og2_genes))
            print(','.join([ogname, og1_count, og2_count]))
