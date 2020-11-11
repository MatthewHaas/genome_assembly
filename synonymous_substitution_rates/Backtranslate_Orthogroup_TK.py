#!/usr/bin/env python
"""Back-translate the protein sequences used to generate the orthologous groups
to nucleotides. This is useful for calculating dN/dS for individual genes in the
orthologue alignment. Requires Biopython and samtools. Takes two arguments:
    1) Directory of FASTA files with nucleotide CDS sequences
    2) Amino acid alignment

There are a couple weird things we have to keep track of for making sure we can
get the sequence from the database. These species have '.p' at the end of their
amino acid sequence names that need to be removed for fetching the CDS:
    - Panicum virgatum
    - Brachypodium distachyon
    - Sorghum bicolor
    - Setaria viridis
    - Setaria_italica
Triticum urartu needs to have '-P1' replaced with '-T1'
B73 and PH207 names end in '_1' and the CDS does not.
"""

import sys
import subprocess
import os
import re # needed for lookup() function
try:
    from Bio import SeqIO
except ImportError:
    print('This script requires Biopython.')
    exit(1)


def list_files(directory):
    """Get the FASTA files that are present in the supplied directory."""
    try:
        abpath = os.path.abspath(directory)
        filepaths = os.listdir(abpath)
        filepaths = [
            os.path.join(abpath, f)
            for f
            in filepaths
            if f.endswith('fa')
            ]
    except (OSError, IOError):
        print('The directory supplied is not readable, or does not exist!')
        exit(1)
    return filepaths

# Lookup definition to go with the LOOKUP table below (part of getting species name).
def lookup(s, lookups):
    for pattern, value in lookups:
        if re.search(pattern, s):
            return value
    return None

def fix_seqname(sname):
    """Perform the 'fixes' for sequence fetching from the indexed CDS files.
    This just makes sure the name matches the ones listed in the FASTA."""
    #   protid is on each line of the FASTA file; splitting doesn't really do anything
    # protid = sname.split(' ')
    # TK 2020-07-22
    # Dictionary for filenames so that we know which CDS file to query for each
    # protein ID.
    lookups = {
        'AET' : 'Aegilops_tauschii.Aet_v4.0.cds.all.fa',
	'PNS' : 'Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa',
	'PNT' : 'Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa',
	'KQJ' : 'Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa',
	'KQK' : 'Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cds.all.fa',
	'Dr' : 'Dioscorea_rotundata.TDr96_F1_Pseudo_Chromosome_v1.0.cds.all.fa',
	'Et' : 'Eragrostis_tef.ASM97063v1.cds.all.fa',
	'HORVU' : 'Hordeum_vulgare.IBSC_v2.cds.all.fa',
	'LPERR' : 'Leersia_perrieri.Lperr_V1.4.cds.all.fa',
	'GSMUA' : 'Musa_acuminata.ASM31385v1.cds.all.fa',
	'OBART' : 'Oryza_barthii.O.barthii_v1.cds.all.fa',
	'ORGLA' : 'Oryza_glaberrima.Oryza_glaberrima_V1.cds.all.fa',
	'ONIVA': 'Oryza_nivara.Oryza_nivara_v1.0.cds.all.fa',
	'ORUFI' : 'Oryza_rufipogon.OR_W1943.cds.all.fa',
	'PVH' : 'Panicum_hallii_fil2.PHallii_v3.1.cds.all.fa',
	'Sspon' : 'Saccharum_spontaneum.Sspon.HiC_chr_asm.cds.all.fa',
	'KQL' : 'Setaria_italica.Setaria_italica_v2.0.cds.all.fa',
	'TraesCS' : 'Triticum_aestivum.IWGSC.cds.all.fa',
	'Zm' : 'Zea_mays.B73_RefGen_v4.cds.all.fa',
	'Zlat': 'Zlat_V1.cds.fa',
        'FUN': 'rice.transcripts.fa',
        'Os': 'Oryza_sativa.IRGSP-1.0.cds.all.fa'
        }
    # Get the filename based on what the sequence starts with.
    for id_start, cds_file in lookups.items():
        if sname.startswith(id_start):
            target_file = cds_file
            break
    # Return the protein name and CDS target file as a tuple
    return (target_file, sname)

    # Make a lookup table to get the species name based on the protein ID.
    # lookups = [('Zlat*','Zizania_latifolia'),('FUN*','Zizania_palustris'),('Os*','Oryza_sativa')]
    # Initialize an empty species dictionary to assist in connecting protid (gene name) to species name
    # species_dict = {}
    # # This for loop will populate the species dictionary so that we can get species name keyed on the protid (gene name)
    # for i in protid:
    #     species = lookup(i, lookups)
    #     return species.encode, i
    #     species_dict[protid] = species.encode()
    # return None


def extract_cds(cds_dir, species, prot_id):
    """Given the protein ID from the alignment and the file list, get the
    correct FASTA to extract from, then use samtools to fetch the CDS sequence.
    """
    #   Get the filename that we want to query
    # for f in file_list:
        #if species in f:
        #    break
    #   Then, buld the command line
    cds_fname = os.path.join(cds_dir, species)
    cmd = ['samtools', 'faidx', cds_fname, prot_id]
    proc = subprocess.Popen(
        cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = proc.communicate()
    #   Then, process the output. We will remove the first line, since it is
    #   just the sequence ID. We will then put all the nucleotides into a single
    #   long string.
    out = out.decode('utf-8')
    lines = out.split('\n')
    cds = ''.join(lines[1:])
    return cds


def backtranslate(p_seq, n_seq):
    """Iterate through the aligned protein sequence, and replace the amino acids
    with codon triplets from the CDS file."""
    #   Keep track of the new sequence. Also keep track of which codon we are
    #   actually processing (gaps don't count)
    newseq = ''
    codon = 0
    for aa in p_seq:
        if aa == '-':
            newseq += '---'
        else:
            newseq += n_seq[codon*3:(codon*3) + 3]
            codon += 1
    return newseq


def main(db_dir, msa):
    """Main function."""
    #   Get the paths of the CDS files
    # paths = list_files(db_dir)
    #   Parse the alignment
    aln = SeqIO.parse(msa, 'fasta')
    for sequence in aln:
        s, p = fix_seqname(sequence.id)
        # TK 2020-07-22
        # Give the path to the database directory to the extract_cds() function
        # because we have already determined which CDS file to query in the
        # fix_seqname() function.
        cdsseq = extract_cds(db_dir, s, p)
        bt_seq = backtranslate(sequence.seq, cdsseq)
        #    Print out the sequence in FASTA format
        print('>' + sequence.id + '\n' + bt_seq)
    return


main(sys.argv[1], sys.argv[2])
