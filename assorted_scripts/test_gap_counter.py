#!/usr/bin/python3.6
import sys

#import the SeqIO module from Biopython
from Bio import SeqIO
with open(sys.argv[1], mode="r") as fasta_handle:
    for record in SeqIO.parse(fasta_handle, "fasta"):
    	start_pos=0
    	counter=0
    	gap=False
    	gap_length = 0
    for char in record.seq:
        if char == 'N':
            if gap_length == 0:
                start_pos=counter
                gap_length = 1
                gap = True
            else:
                gap_length += 1
        else:
            if gap:
                print(record.id + "\t" + str(start_pos) + "\t" + str(start_pos + gap_length))
                gap_length = 0
                gap = False
        counter += 1
