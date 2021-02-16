#!/usr/bin/env python
"""Extract nucleotide CDS from a reference FASTA and GFF file. We will have to
keep track of some weirdness with the phase of the GFF - see the body of the
script for details. Requires Biopython. Takes three arguments:
    1) GFF3
    2) FASTA
    3) Funannotate-GeneName translation table"""

import sys
import re
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation, CompoundLocation
except ImportError:
    sys.stderr.write('This script requires Biopython.\n')
    sys.exit(1)

# Set arguments
try:
    gff_in = sys.argv[1]
    fasta_in = sys.argv[2]
    transt = sys.argv[3]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


def parse_gff(g):
    """Parse the GFF file, and return two dictionaries that describe the
    features and the parent/child relationships."""
    # We also want to store the mRNA->gene information!
    mrna_par = {}
    # And the CDS->mRNA information
    cds_dat = {}
    with open(g, 'r') as f:
        for line in f:
            # if the line is empty or starts with a #, we will skip it
            if line.startswith('#') or line == '\n':
                continue
            else:
                tmp = line.strip().split('\t')
                feat_type = tmp[2]
                if feat_type == 'mRNA':
                    meta = tmp[8].split(';')
                    for m in meta:
                        if m.startswith('ID='):
                            tx_id = m.split('=')[1]
                        if m.startswith('Parent='):
                            tx_par = m.split('=')[1]
                    mrna_par[tx_id] = tx_par
                elif feat_type == 'CDS':
                    scaf = tmp[0]
                    start = tmp[3]
                    end = tmp[4]
                    strand = tmp[6]
                    phase = tmp[7]
                    meta = tmp[8].split(';')
                    for m in meta:
                        if m.startswith('ID='):
                            cds_id = m.split('=')[1]
                        if m.startswith('Parent='):
                            cds_par = m.split('=')[1]
                    if strand == '-':
                        strand = -1
                    else:
                        strand = 1
                    # Watch out for transcripts where there are multiple CDS.
                    # This will require a nested dictionary of lists.
                    if cds_par in cds_dat:
                        pass
                    else:
                        cds_dat[cds_par] = {}
                    if cds_id in cds_dat[cds_par]:
                        pass
                    else:
                        cds_dat[cds_par][cds_id] = []
                    # We want to make a SequenceFeature for each CDS chunk
                    # Keep in mind that GFF is 1-based, so we have to adjust
                    # the start position!
                    cds_feat = SeqFeature(
                        FeatureLocation(int(start)-1, int(end), strand=strand),
                        type="CDS",
                        id=cds_id)
                    # Add some qualifiers to modify the behavior
                    # Use the "standard" genetic code from NCBI
                    cds_feat.qualifiers['transl_tabl'] = [1]
                    # Then, append it into the corresponding dictionary item
                    # keeping the chromosome (scaffold) name and phase with it
                    cds_dat[cds_par][cds_id].append((cds_feat, scaf, phase))
                else:
                    continue
    return (mrna_par, cds_dat)


def combine_features(c_dat):
    """Step through the CDS dictionary and process it down by combining pieces
    that had to be listed as separate chunks in the GFF."""
    # They are keyed on transcript ID
    for tx in c_dat:
        for cds in c_dat[tx]:
            cds_pieces = c_dat[tx][cds]
            # If there fewer than 2 CDS chunks, then pull the tuple out of the
            # list.
            if len(cds_pieces) < 2:
                c_dat[tx][cds] = cds_pieces[0]
            else:
                # Join pieces
                locs = []
                ph = []
                for chunk in cds_pieces:
                    c_loc = FeatureLocation(
                        chunk[0].location.start,
                        chunk[0].location.end,
                        strand=chunk[0].strand)
                    locs.append(c_loc)
                    ph.append(chunk[2])
                # Sort them, according to strand. We assume that a CDS is not a
                # mixed-strand feature
                if cds_pieces[0][0].strand == 1:
                    locs.sort(key=lambda x: x.start)
                else:
                    locs.sort(key=lambda x: x.end, reverse=True)
                # Join them into a CompoundLocation
                full_loc = CompoundLocation(locs)
                # And then overwrite the input dictionary values
                full_feat = SeqFeature(full_loc, type='CDS',
                                       id=cds_pieces[0][0].id)
                full_feat.qualifiers['transl_tabl'] = [1]
                # Keep the phases!
                c_dat[tx][cds] = (full_feat, cds_pieces[0][1], ph)
    return c_dat


def parse_translation(transl):
    """Simple script to return a dictionary keyed on Funannotate ID that stores
    the chosen gene name from the Kimball group."""
    t_table = {}
    with open(transl, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            fun_id = tmp[2]
            gene_name = tmp[0]
            t_table[fun_id] = gene_name
    return t_table


def process_cds(cfs, ref):
    """Extract the coding sequence from the given reference sequence and return
    it."""
    # unpack the tuple
    feat, scaffold, phase = cfs
    # First, extract the sequence of the CDS from the scaffold. This should
    # respect the strand, so we won't have to reverse-complement
    featseq = feat.extract(ref[scaffold])
    return featseq


def main(gff, fa, tt):
    """Main function."""
    # First, we want to parse the GFF3 and store it in memory.
    tx_gene, cds_features = parse_gff(gff)
    # Read in the funannotate-gene name translation table
    fun_gene_table = parse_translation(tt)
    # Then we process the CDS dictionary to make complete CDS features for
    # every transcript
    cds_features = combine_features(cds_features)
    # Read the reference genome into memory
    ref = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    # For each CDS feature, calculate its flags and then print it out. The
    # data are keyed on transcript ID
    cds_seqs = []
    for tx in cds_features:
        for cds in cds_features[tx]:
            cds_nuc = process_cds(cds_features[tx][cds], ref)
            # Get the new gene name from the table. The CDS and mRNA IDs are
            # derived by the new name, so we only have to look up the gene name
            fun_id = tx_gene.get(tx)
            gene_name = fun_gene_table.get(fun_id)
            cds_nuc.id = cds.replace(fun_id, gene_name)
            mrna_id = tx.replace(fun_id, gene_name)
            cds_nuc.description = 'Parent=' + mrna_id + ';Gene=' + gene_name
            cds_seqs.append(cds_nuc)
    # Write it out standard out
    SeqIO.write(cds_seqs, sys.stdout, 'fasta')
    return


main(gff_in, fasta_in, transt)
