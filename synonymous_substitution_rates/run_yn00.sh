#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=15g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e run_yn00.err
#PBS -o run_yn00.out
#PBS -N run_yn00

cd /home/jkimball/haasx092/duplicated_orthogroups

module load paml

for orthogroup in $(cat first_orthogroup_block.txt);
do
cd $orthogroup
yn00 yn00.ctl
cd ..
done
