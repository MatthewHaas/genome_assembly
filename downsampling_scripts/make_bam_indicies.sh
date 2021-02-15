#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=15g,walltime=2:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e make_bam_indicies.err
#PBS -o make_bam_indicies.out
#PBS -N make_bam_indicies

# Purpose of this code is to index downsampled bam files prior to running scythe_mpileup script on the sorted downsampled bam files

cd /home/jkimball/haasx092/pilot_GBS

module load samtools

for file in $(cat downsampled_2fold_bam_list); do
samtools index -b $file
done

for file in $(cat downsampled_4fold_bam_list); do
samtools index -b $file
done

for file in $(cat downsampled_8fold_bam_list); do
samtools index -b $file
done
