#!/bin/bash -l
#PBS -l nodes=1:ppn=24,mem=30g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e run_filter_snps.err
#PBS -o run_filter_snps.out
#PBS -N run_filter_snps

cd /home/jkimball/haasx092/main_GBS/200429_samtools

module load R/3.6.0

Rscript filter_snps_and_make_wide_format.R
