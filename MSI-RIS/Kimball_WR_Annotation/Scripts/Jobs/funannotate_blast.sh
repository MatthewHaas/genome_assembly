#!/bin/bash -l
#PBS -l nodes=1:ppn=24,pmem=1950mb,walltime=36:00:00
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -M mmacchie@umn.edu
#PBS -m abe
#PBS -q mangi

module load ncbi_blast+/2.8.1

# Split the longest protein isoforms into multiple files for blast
# cd /home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/pasa/
# split -l 4000 -a 2 -d augustus.pasa.proteins.longestiso.fa WR_Prot_

cd /home/jkimball/shared/WR_Annotation

BLAST_DB="/panfs/roc/risdb_new/blast/current/nr"
#PROT="/home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/genome.proteins.fa"
# Find all protein FASTA files - uses transcripts IDs
#PROTS=($(find /home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/pasa/PASA_proteins_BLAST/ -mindepth 1 -maxdepth 1 -name 'WR_Prot_*' | sort -V))

# Find all protein FASTA files - uses gene IDs
PROTS=($(find /home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/pasa/PASA_proteins_BLAST2/ -mindepth 1 -maxdepth 1 -name 'WR_Prot_*' | sort -V))
CURR_PROT=${PROTS[${PBS_ARRAYID}]}

# This blast database is from Jan 2019
blastp \
    -query "${CURR_PROT}" \
    -db "${BLAST_DB}" \
    -out "WR_Protein_NR_BLASTP_${PBS_ARRAYID}.xml" \
    -evalue "1e-10" \
    -outfmt 5 \
    -max_target_seqs 3 \
    -qcov_hsp_perc 50 \
    -num_threads "${PBS_NUM_PPN}"


