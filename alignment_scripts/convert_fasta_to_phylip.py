from Bio import SeqIO
import sys

def convert_format(input_file, output_file):
    records = SeqIO.parse(input_file, "fasta")
    count = SeqIO.write(records, output_file, "phylip")
    print("Converted %i records" % count)
    return

convert_format(sys.argv[1], sys.argv[2])
