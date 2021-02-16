#!/bin/bash -l
#PBS -l nodes=1:ppn=24,pmem=1950mb,walltime=6:00:00
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q mangi

module load ncbi_blast+/2.8.1

DB="/panfs/roc/scratch/konox006/JKimball/RNASeq/Uniprot/uniprot_sprot.fasta"
PROT="/home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/predict_results/rice.proteins.fa"

cd /home/jkimball/shared/WR_Annotation/Transcriptome_QC

blastp \
    -query "${PROT}" \
    -db "${DB}" \
    -out Augustus_WR3_Uniprot_blastp.outfmt5 \
    -evalue 1e-20 \
    -num_threads 24 \
    -max_target_seqs 1 \
    -outfmt 5
