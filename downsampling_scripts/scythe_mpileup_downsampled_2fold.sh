#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e scythe_mpileup_2fold.err
#PBS -o scythe_mpileup_2fold.out
#PBS -N scythe_mpileup_2fold

cd /home/jkimball/haasx092/pilot_GBS

module load samtools
module load bcftools
module load htslib
module load parallel

export bams="downsampled_2fold_bam_list"
export prefix="200519_downsampled_2fold"
export ref="/home/jkimball/mshao/genome_seq/zizania_palustris_13Nov2018_okGsv.fasta"

parallel_samtools_processes=15

mkdir -p $prefix
cut -f 1 ${ref}.fai


scythe_mpileup() {
        REGIONS=${1}
        SCAFFOLD=$(echo ${REGIONS} | cut -f 1 -d ";")
        samtools mpileup -q 20 -gDVu \
