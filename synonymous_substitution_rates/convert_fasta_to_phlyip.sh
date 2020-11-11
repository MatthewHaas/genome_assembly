#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e convert_fasta_to_phylip.err
#PBS -o convert_fasta_to_phylip.out
#PBS -N convert_fasta_to_phylip

cd /home/jkimball/haasx092/duplicated_orthogroups

module load python

for orthogroup in $(cat first_orthogroup_block.txt); 
do
python convert_fasta_to_phylip.py $orthogroup/${orthogroup}_slim_backtranslated.fa $orthogroup/${orthogroup}_slim_backtranslated.phy
done
