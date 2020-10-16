#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=50g,walltime=8:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e clustal_omega.err
#PBS -o clustal_omega.out
#PBS -N clustal_omega

# This script is meant to make an aligned peptide file for 4DTV analysis

cd /home/jkimball/haasx092/4DTV_latifolia

module load clustalo/1.2.1_gcc4.9.2_haswell

clustalo -i OG0000003.fa  -o clustal_omega_output.fa
