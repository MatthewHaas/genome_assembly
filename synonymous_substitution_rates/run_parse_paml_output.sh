#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=15g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e parse_paml.err
#PBS -o all_yn00_data.txt
#PBS -N parse_paml

cd /home/jkimball/haasx092/duplicated_orthogroups/zizania_specific_duplications

module load python3

for orthogroup in $(cat zizania_duplications.txt);
do
python parse_paml_output.py $orthogroup/yn00_${orthogroup}_output_file.txt
done
