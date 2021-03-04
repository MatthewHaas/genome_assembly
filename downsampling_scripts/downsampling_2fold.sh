#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=10:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o downsampling.out
#SBATCH -e downsampling.err

# This code is for downsampling fastq files (2-fold)

cd /home/jkimball/haasx092/pilot_GBS

for i in $(cat fastq_list); do
STEM=$(echo ${i} | cut -f 1 -d "_")
 zcat ${i} | awk '!(int((NR - 1)/4) % 2)' | gzip > ${STEM}_2fold.fq.gz 
done
