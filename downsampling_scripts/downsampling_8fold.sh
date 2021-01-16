#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e downsampling_8fold.err
#PBS -o downsampling_8fold.out
#PBS -N downsampling_8fold

# This code is for downsampling fastq files (8-fold)

cd /home/jkimball/haasx092/pilot_GBS

for i in $(cat fastq_list); do
STEM=$(echo ${i} | cut -f 1 -d "_")
 zcat ${i} | awk '!(int((NR - 1)/4) % 8)' | gzip > ${STEM}_8fold.fq.gz 
done
