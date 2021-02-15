#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=15g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e bwa_2fold_down.err
#PBS -o bwa_2fold_down.out
#PBS -N bwa_2fold_down

# Purpose of this script is to align downsampled (2-fold) fastq files to the NWR genome for SNP calling

cd /home/jkimball/haasx092/pilot_GBS

module load bwa
module load samtools

FASTA='/home/jkimball/mshao/genome_seq/zizania_palustris_13Nov2018_okGsv.fasta.gz'

for i in $(cat downsampling_2fold.txt); do
STEM=$(echo ${i} | cut -f 1 -d "/")
bwa mem $FASTA ${i} 2> $STEM/${STEM}_downsampling_2fold_bwa.err | samtools sort -o $STEM/${STEM}_downsampling_2fold_sorted.bam 2> $STEM/${STEM}_downsampling_2fold_samtools_sort.err;
done
