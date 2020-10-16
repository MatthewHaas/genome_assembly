#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=15g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e normalize_awk.err
#PBS -o normalize_awk.out
#PBS -N normalize_awk

cd /home/jkimball/haasx092/pilot_GBS/200511_samtools

# This was done inside the samtools directory
normalize_prefix="190511_normalize"
zcat ${prefix}_*.vcf.gz | awk -f ./normalize.awk > $normalize_prefix.tsv 2> $normalize_prefix.err
