#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e LTR_density_calc.err
#PBS -o LTR_density_calc.out
#PBS -N LTR_density_calc

cd /home/jkimball/haasx092/circos/200511_version

module load R/3.6.0

Rscript make_LTR_density_files_for_circos.R