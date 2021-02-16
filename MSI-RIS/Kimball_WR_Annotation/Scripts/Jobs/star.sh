#!/bin/bash
#PBS -l mem=250GB,nodes=1:ppn=24,walltime=35:00:00
#PBS -A jkimball 
#PBS -m abe
#PBS -j oe
#PBS -M mmacchie@umn.edu
#PBS -q ram256g
#PBS -W group_list=jkimball
#PBS -N star

module load star/2.7.1a



#GENOME='/panfs/roc/scratch/konox006/JKimball/RNASeq/genome/'
#R1='/panfs/roc/scratch/konox006/JKimball/RNASeq/Trimmed_reads/WR_8Tissues.R1.trimmed.fastq'
#R2='/panfs/roc/scratch/konox006/JKimball/RNASeq/Trimmed_reads/WR_8Tissues.R2.trimmed.fastq'
#INDEX='/panfs/roc/scratch/konox006/JKimball/RNASeq/genome/'
#OUT='/panfs/roc/scratch/konox006/JKimball/RNASeq/Aligned_reads/'

echo -n "Ran: "
date


# Index rice genome
#STAR --runMode genomeGenerate --runThreadN 16 --genomeFastaFiles ${GENOME} --genomeDir ${INDEX}

# did not work
#STAR --runMode alignReads --runThreadN 24 --genomeDir ${INDEX} --readFilesIn ${R1} ${R2} --outFileNamePrefix ${OUT}/rice_rnas --outSAMtype BAM Unsorted

# Align all reads - works
#STAR --runMode alignReads --runThreadN 24 --genomeDir ${INDEX} --readFilesIn ${R1} ${R2} --outFileNamePrefix ${OUT}/rice_rnas2 --limitBAMsortRAM 24846526061 --outSAMtype BAM SortedByCoordinate


GENOME='/home/jkimball/shared/WR_Annotation/genome/reformatted_genome/zizania_palustris_13Nov2018_okGsv.fasta.masked.headerreformatted'
R1='/home/jkimball/shared/WR_Annotation/Trimmed_reads/WR_8Tissues.R1.trimmed.fastq'
R2='/home/jkimball/shared/WR_Annotation/Trimmed_reads/WR_8Tissues.R2.trimmed.fastq'
OUT='/home/jkimball/shared/WR_Annotation/Aligned_reads/'
INDEX='/home/jkimball/shared/WR_Annotation/genome/reformatted_genome/'
STAR --runMode genomeGenerate --runThreadN 24 --genomeFastaFiles ${GENOME} --genomeDir ${INDEX}
STAR --runMode alignReads --runThreadN 24 --genomeDir ${INDEX} --readFilesIn ${R1} ${R2} --outFileNamePrefix ${OUT}/rice_rnas_chromreformatted --limitBAMsortRAM 24846526061 --outSAMtype BAM SortedByCoordinate


echo -n "Done: "
date



