#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e convert_fasta_to_phylip.err
#PBS -o convert_fasta_to_phylip.out
#PBS -N convert_fasta_to_phylip

cd /home/jkimball/haasx092/Ks_substitution_rates

module load python

python convert_fasta_to_phylip.py test_orthopair_aligned.fa output_file.nuc
