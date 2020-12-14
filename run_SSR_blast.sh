#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=15g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e blast_SSR.err
#PBS -o blast_SSR.out
#PBS -N blast_SSR

# 13 December 2020
# Make sure you're in the right directory...
cd /home/jkimball/haasx092/align_SSR_markers

# Load BLAST
module load ncbi_blast+

# This loop iterates over each line in the regular text file "shattering_fasta_list.txt" and performs a BLAST each time.
# The BLAST is to look for orthologs of shattering genes in other plant species (with an emphasis on rice/grasses)
# Each line represents a different FASTA file.

for i in $(cat marker_list.txt)
do
STEM=$(echo ${i} | cut -f 1 -d ".")
blastn -num_threads 10 -db ~/genome_seq/zizania_palustris_genome_db -query $i > ${STEM}_out.txt
done
