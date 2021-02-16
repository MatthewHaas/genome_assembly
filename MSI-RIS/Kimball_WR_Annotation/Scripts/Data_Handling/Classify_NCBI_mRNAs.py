#!/usr/bin/env python
"""Classify the mRNA features in the Northern Wild Rice genome annotation
according to whether they have identical protein products and different UTRs.
Takes one argument:
    1) Feature Table file to process (gzipped)"""

import sys
import gzip

try:
    ft_gz = sys.argv[1]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)


def parse_feat(ft):
    """Read through the feature table and return a dictionary of the genes and
    mRNAs. We should also keep the positions so we can check if the UTRs are
    identical and/or the CDS is identical."""
    ft_d = {}
    # We have to open the file in this strange way because we need to control
    # how we iterate through the file in a way that is incompatible with a
    # regular for loop.
    ft_handle = gzip.open(ft, 'rt')
    ft_iter = iter(ft_handle)
    for line in ft_iter:
        # If the entry ends in "gene" then we have found a gene for us to
        # process
        if line.endswith('gene\n'):
            start, stop, dummy = line.strip().split('\t')
            # Advance the iterator to get the locus tag, too
            dummy, loc_tag = next(ft_iter).strip().split('\t')
            # Push the gene into the dictionary
            ft_d[loc_tag] = {}
        # If the entry ends in "mRNA" then we found a transcript to process
        if line.endswith('mRNA\n'):
            starts = []
            stops = []
            mrna_line = line
            while 'transcript_id' not in mrna_line:
                if 'mRNA' in mrna_line:
                    start, stop, dummy = mrna_line.strip().split('\t')
                else:
                    start, stop = mrna_line.strip().split('\t')
                starts.append(start)
                stops.append(stop)
                mrna_line = next(ft_iter)
            if 'transcript_id' in mrna_line:
                txid = mrna_line.strip().split('\t')[1]
                parent_id = txid.split('|')[-1]
                parent_id = parent_id.split('-')[0]
                # Make intervals out of the starts and stops
                intervals = tuple(zip(starts, stops))
                # Push this data into the dictionary
                ft_d[parent_id][txid] = {'intervals': intervals}
        # We will also be very trusting with the CDS annotations
        if line.endswith('CDS\n'):
            starts = []
            stops = []
            cds_line = line
            while 'transcript_id' not in cds_line:
                if 'CDS' in cds_line:
                    start, stop, dummy = cds_line.strip().split('\t')
                else:
                    start, stop = cds_line.strip().split('\t')
                starts.append(start)
                stops.append(stop)
                cds_line = next(ft_iter)
            # Then, we have to process these a bit differently from the mRNA.
            # We will take the next indented blocks as the meta lines for the
            # CDS.
            if 'transcript_id' in cds_line:
                txid = cds_line.strip().split('\t')[1]
                gene_id = txid.split('|')[-1].split('-')[0]
                intervals = tuple(zip(starts, stops))
                # Push this into the dictionary, too
                ft_d[gene_id][txid]['CDS'] = intervals
    ft_handle.close()
    return ft_d


def process_gene(ft_gene):
    """Process the parsed gene feature with its mRNA and CDS children features.
    """
    # If there is only one transcript, then we don't need to adjust anything!
    if len(ft_gene) == 1:
        return None
    # If there are multiple transcripts, we will then compare their coding
    # sequences.
    identical_cds = {}
    for mrna in ft_gene:
        if ft_gene[mrna]['CDS'] in identical_cds:
            identical_cds[ft_gene[mrna]['CDS']].append(mrna)
        else:
            identical_cds[ft_gene[mrna]['CDS']] = [mrna]
    # Then, within each group of transcripts that have the same CDS, we want to
    # compare the UTRs
    identical_utrs = {}
    for prod in identical_cds:
        identical_utrs[prod] = {}
        for mrna in identical_cds[prod]:
            utr = []
            for i in ft_gene[mrna]['intervals']:
                if i in prod:
                    continue
                else:
                    utr.append(i)
            # cast to tuple so we can use it as a dictionary key
            utr = tuple(utr)
            # And push it into the dictionary
            if utr in identical_utrs[prod]:
                identical_utrs[prod][utr].append(mrna)
            else:
                identical_utrs[prod][utr] = [mrna]
    # OK, yikes, we have a big ugly data structure now. Let's re-organize it
    # for easy printing, sorting on CDS length.

    def cdslen(intervals):
        """Return the length of a feature calculated from a list of tuples of
        (start, stop) values"""
        feat_len = 0
        for i in intervals:
            sane_start = i[0].translate({ord(x): None for x in '<>'})
            sane_end = i[1].translate({ord(x): None for x in '<>'})
            sub_len = int(sane_end) - int(sane_start)
            feat_len += sub_len
        return feat_len

    organized = {}
    prod_type = 1
    for cds in sorted(identical_utrs, key=cdslen):
        new_key = 'Type_' + str(prod_type)
        organized[new_key] = {}
        tx_variant = 1
        for utr in sorted(identical_utrs[cds]):
            new_var = 'Variant_' + str(tx_variant)
            if new_var in organized[new_key]:
                organized[new_key][new_var] += identical_utrs[cds][utr]
            else:
                organized[new_key][new_var] = identical_utrs[cds][utr]
            tx_variant += 1
        prod_type += 1
    return organized


def main(feat_table):
    """Main function."""
    # First we want to read through the feature table and store the information
    # about the genes and mRNAs
    ft_dat = parse_feat(feat_table)
    # Now we want to process the parsed feature table
    # Print a header!
    print(','.join(['Gene.ID', 'Transcript.ID', 'CDS.Type', 'UTR.Type']))
    # Also keep a dictionary for the slicing state of the transcripts
    trans_splice = {}
    for gene in sorted(ft_dat):
        gene_info = process_gene(ft_dat[gene])
        if not gene_info:
            # Ugly way to get the transcript ID
            txid = list(ft_dat[gene].keys())[0]
            toprint = [
                gene,
                txid,
                '1',
                '1']
            print(','.join(toprint))
            trans_splice[txid] = 'NoAlt'
        else:
            for cds_type in sorted(gene_info, key=lambda x: int(x.split('_')[1])):
                # Check the length of the list of transcripts with the same
                # CDS. If it is greater than 1, then the same protein
                # product is coded for by different transcripts.
                if len(gene_info[cds_type]) > 1:
                    splice = 'SameProduct'
                else:
                    splice = 'DifferentProduct'
                for utr_type in sorted(gene_info[cds_type], key=lambda x: int(x.split('_')[1])):
                    txid = gene_info[cds_type][utr_type][0]
                    toprint = [
                        gene,
                        txid,
                        cds_type.split('_')[1],
                        utr_type.split('_')[1]]
                    print(','.join(toprint))
                    trans_splice[txid] = splice
    # Now, print the splicing data to stderr
    sys.stderr.write('Transcript.ID,Splicing\n')
    for transcript in sorted(trans_splice):
        towrite = [transcript, trans_splice[transcript]]
        sys.stderr.write(','.join(towrite) + '\n')
    return


main(ft_gz)
