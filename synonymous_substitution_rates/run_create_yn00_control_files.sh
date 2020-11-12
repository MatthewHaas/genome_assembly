#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e make_yn00_control_files.err
#PBS -o make_yn00_control_files.out
#PBS -N make_yn00_control_files

cd /home/jkimball/haasx092/duplicated_orthogroups

module load python3

for orthogroup in $(cat first_orthogroup_block.txt);
do
python create_yn00_control_files.py yn00_template.ctl $orthogroup/yn00.ctl $orthogroup
done
