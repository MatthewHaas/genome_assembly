#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=4:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o index_bams.out
#SBATCH -e index_bams.err

# Purpose of this code is to index downsampled bam files prior to running scythe_mpileup script on the sorted downsampled bam files

cd /home/jkimball/haasx092/pilot_GBS

module load samtools

for file in $(cat downsampled_2fold_bam_list); do
samtools index -b $file
done

for file in $(cat downsampled_4fold_bam_list); do
samtools index -b $file
done

for file in $(cat downsampled_8fold_bam_list); do
samtools index -b $file
done
