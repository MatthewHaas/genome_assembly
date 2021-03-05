#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=30g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e backtranslate.err
#PBS -o backtranslate.out
#PBS -N backtranslate

# First go to the directory with aligned sequences to make the file with list of sequences
cd /home/jkimball/haasx092/duplicated_orthogroups

module load samtools
module load python3

for orthogroup in $(cat first_orthogroup_block.txt); do
python Backtranslate_Orthogroup_TK.py cds_fasta_files $orthogroup/${orthogroup}_slim_aligned.fa > $orthogroup/${orthogroup}_slim_backtranslated.fa
done
