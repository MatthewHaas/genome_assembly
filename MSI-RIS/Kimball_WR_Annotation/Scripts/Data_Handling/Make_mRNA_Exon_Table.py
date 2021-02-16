#!/usr/bin/env python
"""Digest a GFF and print out information about mRNA and exons for a NCBI
feature table. This will account for the strand of the transcripts and their
parent gene features. Takes one argument:
    1) GFF"""

import sys

try:
    gff_in = sys.argv[1]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


def parse_gff(g):
    """Parse a GFF and return a dictionary that stores the transcript
    relationships with its parent gene and its child exons. We really only
    care about the positions and the strand."""
    gff_dat = {}
    with open(g, 'r') as f:
        for line in f:
            # Ignore comment/empty lines
            if line.startswith('#') or line == '\n':
                continue
            else:
                tmp = line.strip().split('\t')
                scaffold = tmp[0]
                # Let's convert the strand into the same convention used by
                # the CDS script, for consistency
                if tmp[6] == '+':
                    strand = '1'
                else:
                    strand = '-1'
                meta = tmp[8]
                # First look for mRNA - these should come before exons in
                # the GFF file, so it should be safe to run this check first
                feat_type = tmp[2]
                if feat_type == 'mRNA':
                    # Extract the transcript start and end
                    tx_start = tmp[3]
                    tx_end = tmp[4]
                    # extract the transcript ID
                    for m in meta.split(';'):
                        if m.startswith('ID='):
                            tx_id = m.split('=')[1]
                        if m.startswith('Parent='):
                            g_par = m.split('=')[1]
                    # Put it into the dictionary we are building
                    if tx_id in gff_dat:
                        continue
                    else:
                        gff_dat[tx_id] = {
                            'Parent': g_par,
                            'Strand': strand,
                            'Scaffold': scaffold,
                            'T_Start': tx_start,
                            'T_End': tx_end}
                elif feat_type == 'exon':
                    # For exons, extract the start and end
                    e_start = tmp[3]
                    e_end = tmp[4]
                    # Figure out which transcript this exon is from
                    for m in meta.split(';'):
                        if m.startswith('Parent='):
                            tx_par = m.split('=')[1]
                    # If we have already seen an exon from this transcript,
                    # append it to the list
                    if 'Exons' in gff_dat[tx_par]:
                        gff_dat[tx_par]['Exons'].append((e_start, e_end))
                    else:
                        gff_dat[tx_par]['Exons'] = [(e_start, e_end)]
                else:
                    continue
    return gff_dat


def main(gff):
    """Main function."""
    gdat = parse_gff(gff)
    # Print a header
    header = [
        'Transcript.ID',
        'Gene.ID',
        'Scaffold',
        'mRNA.Start',
        'mRNA.End',
        'Exon.Starts',
        'Exon.Ends',
        'Strand']
    print('\t'.join(header))
    # Then, print out the table
    for tx in sorted(gdat):
        gpar = gdat[tx]['Parent']
        scaf = gdat[tx]['Scaffold']
        tx_start = gdat[tx]['T_Start']
        tx_end = gdat[tx]['T_End']
        strand = gdat[tx]['Strand']
        # Rearrange the start/end based on the strand
        if strand == '1':
            ft_mrna_start = tx_start
            ft_mrna_end = tx_end
        else:
            ft_mrna_start = tx_end
            ft_mrna_end = tx_start
        ft_exon_starts = [
            x[0]
            if strand == '1'
            else x[1]
            for x
            in gdat[tx]['Exons']]
        ft_exon_ends = [
            x[1]
            if strand == '1'
            else x[0]
            for x
            in gdat[tx]['Exons']]
        # If we are looking at a reversely-stranded gene, then we have re-sort
        # the positions in descending order
        if strand == '-1':
            ft_exon_starts.sort(key=lambda x: int(x), reverse=True)
            ft_exon_ends.sort(key=lambda x: int(x), reverse=True)
        # Join start and end positions on a comma for printing
        p_starts = ','.join(ft_exon_starts)
        p_ends = ','.join(ft_exon_ends)
        # And then print it out
        toprint = [tx, gpar, scaf, ft_mrna_start, ft_mrna_end, p_starts,
                   p_ends, strand]
        print('\t'.join(toprint))
    return


main(gff_in)
