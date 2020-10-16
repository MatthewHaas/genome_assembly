#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=15g,walltime=1:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e test_gap_counter.err
#PBS -o test_gap_counter.out
#PBS -N test_gap_counter

cd /home/jkimball/haasx092/genome_seq

module load python3

python test_gap_counter.py Scaffold_1.fa
