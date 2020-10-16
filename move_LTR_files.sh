#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=15g,walltime=0:20:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e move_LTR_files.err
#PBS -o move_LTR_files.out
#PBS -N move_LTR_files

# Purpose of this code is to move various long terminal repeat (LTR) density files to subdirectories to help make circos plot for repeats.
# The files moved by this script come from the shell script run_make_LTR_density_files.sh which runs make_LTR_density_files_for_circos.R.

cd /home/jkimball/haasx092/circos/200511_version

mkdir Copia
mkdir Gypsy
mkdir otherLTRs

for i in $(cat largest_scaffolds.txt); do
mv ./${i}_Copia_binned.txt ./Copia/${i}_Copia_binned.txt
mv ./${i}_Gypsy_binned.txt ./Gypsy/${i}_Gypsy_binned.txt
mv ./${i}_otherLTR_binned.txt ./otherLTRs/${i}_otherLTRs_binned.txt
done