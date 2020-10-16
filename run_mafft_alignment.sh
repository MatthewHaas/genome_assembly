#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e run_mafft_align.err
#PBS -o run_mafft_align.out
#PBS -N mafft_align

cd /home/jkimball/haasx092/4DTV_pascal

#mkdir Aligned_Orthogroup_Seqs

module load mafft
mafft --keeplength Orthogroup_Sequences/OG0000001_trimmed.fa > OG0000001_aligned.fa
#for i in $(cat Orthogroup_Sequences/orthogroup_sequence_list.txt)
#do
#	STEM=$(echo ${i} | cut -f 1 -d ".")
#	mafft Orthogroup_Sequences/${STEM}_trimmed.fa > Aligned_Orthogroup_Seqs/${STEM}_aligned.fa
#done
