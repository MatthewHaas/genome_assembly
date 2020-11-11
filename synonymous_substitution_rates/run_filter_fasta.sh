#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e filter_fasta.err
#PBS -o filter_fasta.out
#PBS -N filter_fasta

cd /home/jkimball/haasx092/duplicated_orthogroups

module load python3

for file in $(cat first_orthogroup_block.txt);
do
python filter_fasta.py $file/${file}.fa $file/${file}_slim.fa
done
