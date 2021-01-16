#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=15g,walltime=0:20:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e make_gap_table.err
#PBS -o make_gap_table.out
#PBS -N make_gap_table

# This script counts the number of gaps (lines) in each file containing gaps found by the gap counter scripts gap_counter.sh/gap_counter.py

cd /home/jkimball/haasx092/genome_seq

awk 'BEGIN{OFS="\t"; print "Scaffold", "Number_of_gaps", "Original_file"}' gap_count_table.txt >> gap_count_table.txt
for file in $(cat scaffold_gap_files.txt);
do
scaffold=$(echo ${file} | cut -f 1,2 -d "_")
echo -e $scaffold'\t'$(wc -l $file) 1>> gap_count_table.txt
done
sed -i -e 's/ /\t/g' gap_count_table.txt
