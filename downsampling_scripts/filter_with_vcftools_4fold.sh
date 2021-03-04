#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=2:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o filter_vcfs.out
#SBATCH -e filter_vcfs.err

# Include path to the working directory
cd ~/pilot_GBS/200519_downsampled_4fold

# Filter VCF files to keep samples with fewer than 200 genotypes missing (to make a reference panel)
for i in $(cat vcf_list_new_4fold.txt);
do
STEM=$(echo ${i} | cut -f 1 -d ".")
~/vcftools/bin/vcftools --gzvcf  $i --max-missing 0.9 --min-alleles 2 --max-alleles 2 --remove-indels --minDP 6 --recode --recode-INFO-all --out ${STEM}_filtered
done
