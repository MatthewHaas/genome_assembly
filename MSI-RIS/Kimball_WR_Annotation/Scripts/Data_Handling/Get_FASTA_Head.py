#!/usr/bin/env python
"""Get the top N sequences of a FASTA file. Simple script for testing. Requires
Biopython, and takes two arguments:
    1) FASTA to subset
    2) Number of sequences"""

try:
    import sys
    from Bio import SeqIO
    fa_in = sys.argv[1]
    nseqs = int(sys.argv[2])
    assert nseqs > 0
except ImportError:
    sys.stderr.write('This script requires Biopython.\n')
    sys.exit(1)
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(2)
except (ValueError, AssertionError) as e:
    sys.stderr.write('Second argument must be a positive integer.\n')
    sys.exit(3)

seq_h = open(fa_in, 'r')
for index, rec in enumerate(SeqIO.parse(seq_h, 'fasta')):
    if index < nseqs:
        SeqIO.write(rec, sys.stdout, 'fasta')
    else:
        break
seq_h.close()
