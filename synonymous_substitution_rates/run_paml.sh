#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=15g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e run_paml.err
#PBS -o run_paml.out
#PBS -N run_paml

cd /home/jkimball/haasx092/Ks_substitution_rates

module load paml

codeml codeml.ctl
