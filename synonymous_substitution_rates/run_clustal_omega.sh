#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=50g,walltime=8:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e clustal_omega.err
#PBS -o clustal_omega.out
#PBS -N clustal_omega

# This script is meant to make an aligned peptide file for 4DTV analysis

cd /home/jkimball/haasx092/duplicated_orthogroups

#module load clustalo/1.2.1_gcc4.9.2_haswell
module load clustalo/1.2.1_gcc4.9.2


for orthogroup in $(cat first_orthogroup_block.txt);
do
clustalo -i $orthogroup/${orthogroup}_slim.fa -t Protein --infmt=fasta --full-iter --outfmt=fasta -o $orthogroup/${orthogroup}_slim_aligned.fa
done
