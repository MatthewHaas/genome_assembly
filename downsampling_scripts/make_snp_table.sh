#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e make_snp_table.err
#PBS -o make_snp_table.out
#PBS -N make_snp_table

cd /home/jkimball/haasx092/pilot_GBS/200519_downsampled_8fold

prefix="200519_downsampled_8fold"
normalize_prefix="200529_normalize"
mktemp | read tmp
zcat ${prefix}_*.vcf.gz | awk -f ./normalize.awk > '$tmp' 2> $normalize_prefix.err

# Copy the results of SNP calling in a temporary file to a TSV file for further processing

cp '$tmp' 200529_normalize_8fold.tsv # Can further process/filter the TSV file in R.

rm '$tmp' # Remove temporary file
