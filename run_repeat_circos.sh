#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=15g,walltime=0:30:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e run_repeat_circos.err
#PBS -o run_repeat_circos.out
#PBS -N run_repeat_circos

# 12 May 2020
# Circos plot for Zizania palustris

# Change into directory
cd /home/jkimball/haasx092/circos/200511_version

# Load required software (version is important!)
module load circos/0.69-6
module load perl/5.26.0

# Run circos by calling configuration file.
circos -conf repeat_specific_circos.conf -outputfile repeat_specific_circos.png
