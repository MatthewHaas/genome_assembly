#!/bin/bash -l
#PBS -l nodes=1:ppn=24,pmem=1950mb,walltime=48:00:00
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q mangi

module load ncbi_blast+/2.8.1

ASM="/home/jkimball/shared/WR_Annotation/Transcriptome_Asm_trinity/Trinity.fasta"
DB="/panfs/roc/scratch/konox006/JKimball/RNASeq/Uniprot/uniprot_sprot.fasta"
PROT="/home/jkimball/shared/WR_Annotation/Annotation/wild_rice2/predict_misc/augustus_proteins.fa"

cd /home/jkimball/shared/WR_Annotation/Transcriptome_QC

#blastx \
#    -query "${ASM}" \
#    -db "${DB}" \
#    -out Trinity_Uniprot_blastx.outfmt6 \
#    -evalue 1e-20 \
#    -num_threads 24 \
#    -max_target_seqs 1 \
#    -outfmt 6\

blastp \
    -query "${PROT}" \
    -db "${DB}" \
    -out Augustus_Uniprot_blastp.outfmt5 \
    -evalue 1e-20 \
    -num_threads 24 \
    -max_target_seqs 1 \
    -outfmt 5
