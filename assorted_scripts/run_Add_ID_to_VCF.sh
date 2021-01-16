#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e Add_IDs_to_VCF_files.err
#PBS -o Add_IDs_to_VCF_files.out
#PBS -N Add_IDs_to_VCF_files

cd /home/jkimball/haasx092/main_GBS/200305_samtools

module load python3

# Call script. First argument is file name, second argument (zp) is the prefix (can change to be anything), third argument is the number of zeros to pad the SNP identifier
for $vcf_filename in $(cat vcf_list.txt); do
STEM=$(echo ${i} | cut -f 1 -d ".")
python Add_ID_to_VCF.py vcf_filename zp 4 > ${STEM}_with_IDs.vcf
done
