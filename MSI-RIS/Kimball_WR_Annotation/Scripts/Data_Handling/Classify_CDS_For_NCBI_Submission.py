#!/usr/bin/env python
"""Make a lookup table of CDS features for each gene in the Northern Wild Rice
genome annotation that contains information for making the NCBI feature table.
This script will extract the sequence and check it for completeness: length is
a multiple of three, first codon is Met, final codon is Stop. Requires
Biopython. Takes two arguments:
    1) GFF3
    2) FASTA"""

import sys
import gzip
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


def process_cds(cfs, ref):
    """Extract the coding sequence from the given reference sequence and
    process it. This is where we will apply our rules for length/etc."""
    # unpack the tuple
    feat, scaffold, phase = cfs
    # First, extract the sequence of the CDS from the scaffold. This should
    # respect the strand, so we won't have to reverse-complement
    featseq = feat.extract(ref[scaffold])
    # Calculate the length of the codon-containing sequence
    cdslen = len(featseq)
    # Run a check to see if it is a multiple of 3 or not
    if cdslen % 3 == 0:
        mult_three = True
    else:
        mult_three = False
    # We want to take the phase of the first piece of the CDS and remove that
    # number of bases from the beginning. For some reason, it looks like the
    # phase is messed up here - 2 should be 1 and 1 should be 2?
    first_phase = int(phase[0])
    if first_phase == 2:
        first_phase = 1
    elif first_phase == 1:
        first_phase = 2
    # Save the updated phases
    new_phases = []
    for idx, p in enumerate(phase):
        if idx == 0:
            new_phases.append(str(first_phase))
        else:
            new_phases.append(p)
    featseq = featseq[first_phase:]
    # Then we will translate it
    cds_trans = featseq.translate(table=1)
    # How many parts does the CDS have?
    cds_numparts = str(len(feat.location.parts))
    # Extract the start positions and end positions
    cds_starts = [str(x.start) for x in feat.location.parts]
    cds_ends = [str(x.end) for x in feat.location.parts]
    # Check for 5-prime partial CDS
    if str(cds_trans.seq).startswith('M'):
        fiveprime_part = False
    else:
        fiveprime_part = True
    # Check for 3-prime partial CDS
    if str(cds_trans.seq).endswith('*'):
        threeprime_part = False
    else:
        threeprime_part = True
    # Let's do a thing where we check if the next triplet would be a stop
    # codon.
    if threeprime_part:
        # Adjust the positions by one codon, according to the strand
        if feat.strand == 1:
            newstart = feat.location.start
            newend = feat.location.end + 3
        else:
            newstart = feat.location.start - 3
            newend = feat.location.end
        # Make a new feature and set the qaulifiers
        new_feat = SeqFeature(
            FeatureLocation(newstart, newend, strand=feat.strand),
            type="CDS",
            id=feat.id)
        new_feat.qualifiers['transl_tabl'] = [1]
        # Translate it
        new_trans = new_feat.translate(ref[scaffold])
        # Check if it ends in a stop
        if str(new_trans.seq).endswith('*'):
            next_codon_is_stop = True
        else:
            next_codon_is_stop = False
    else:
        next_codon_is_stop = 'NA'
    # Lastly, check for internal stop codons
    internal_stop_re = re.compile(r'^[^\*].+\*.+[^\*]$')
    if internal_stop_re.match(str(cds_trans.seq)):
        has_i_stop = True
    else:
        has_i_stop = False
    cds_flags = (mult_three, fiveprime_part, threeprime_part, has_i_stop,
                 cds_trans, cds_numparts, next_codon_is_stop, cdslen,
                 cds_starts, cds_ends, new_phases)
    return cds_flags


def main(gff, fa):
    """Main function."""
    # First, we want to parse the GFF3 and store it in memory.
    tx_gene, cds_features = parse_gff(gff)
    # Then we process the CDS dictionary to make complete CDS features for
    # every transcript
    cds_features = combine_features(cds_features)
    # Read the reference genome into memory
    ref = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    # Print a header!
    header = ['CDS.ID', 'Chunk.Num.Parts', 'CDS.Len', 'Tx.ID', 'Gene.ID',
              'Scaffold', 'Starts', 'Ends', 'Strand', 'Phase', 'Mult.Three',
              '5prime.Partial', '3prime.Partial', 'Next.Codon.is.Stop',
              'Has.Interal.Stops', 'CDS.Trans']
    print('\t'.join(header))
    # For each CDS feature, calculate its flags and then print it out. The
    # data are keyed on transcript ID
    for tx in cds_features:
        for cds in cds_features[tx]:
            # First, get the parent info
            parent_gene = tx_gene.get(tx, 'Unknown')
            # Classify the sequences
            cds_class = process_cds(cds_features[tx][cds], ref)
            phs = ','.join(cds_class[10])
            # Then print out the info! We will print it out in this order:
            #   - CDS ID
            #   - Number of CDS parts
            #   - CDS length
            #   - Parent transcript ID
            #   - Parent gene ID
            #   - Scaffold
            #   - Start
            #   - End
            #   - Strand
            #   - Multiple of 3?
            #   - 5-prime partial?
            #   - 3-prime partial?
            #   - Next codon is stop?
            #   - Internal stops?
            #   - Translated sequence
            cds_id = cds_features[tx][cds][0].id
            num_parts = cds_class[5]
            cds_len = str(cds_class[7])
            scaffold = cds_features[tx][cds][1]
            strand = str(cds_features[tx][cds][0].strand)
            transeq = str(cds_class[4].seq)
            internal_stop = str(cds_class[3])
            threeprime_part = cds_class[2]
            fiveprime_part = cds_class[1]
            mult_three = str(cds_class[0])
            next_is_stop = str(cds_class[6])
            # This is going to be a bit hairy! We will have to set the "real"
            # start and end based on the GFF positions and the reported strand
            # Then, if a feature is partial in one direction, we will have to
            # modify its start/end, but also keep in mind the strand, i.e., for
            # three-prime partial features, the *first* (sorted by coordinate)
            # feature will have to be modified. Also remember to add 1 to the
            # gff_starts values to convert it back to 1-based coordinates!
            gff_starts = [str(int(i)+1) for i in cds_class[8]]
            gff_ends = cds_class[9]
            if strand == '1':
                real_starts = gff_starts
                real_ends = gff_ends
            else:
                real_starts = gff_ends
                real_ends = gff_starts

            # Check five-prime partialness.
            if fiveprime_part:
                real_starts[0] = '<' + real_starts[0]
            # And the threeprime partialness
            if threeprime_part:
                real_ends[-1] = '>' + real_ends[-1]
            # And turn the starts and stops into printable things
            p_starts = ','.join(real_starts)
            p_ends = ','.join(real_ends)
            toprint = [
                cds_id,
                num_parts,
                cds_len,
                tx,
                parent_gene,
                scaffold,
                p_starts,
                p_ends,
                strand,
                phs,
                mult_three,
                str(fiveprime_part),
                str(threeprime_part),
                next_is_stop,
                internal_stop,
                transeq]
            print('\t'.join(toprint))
    return


main(gff_in, fasta_in)
