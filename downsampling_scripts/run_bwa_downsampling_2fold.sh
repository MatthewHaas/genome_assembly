#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=10:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o run_bwa_2fold.out
#SBATCH -e run_bwa_2fold.err

# Purpose of this script is to align downsampled (2-fold) fastq files to the NWR genome for SNP calling

cd /home/jkimball/haasx092/pilot_GBS

module load bwa
module load samtools

FASTA='/home/jkimball/shared/WR_Annotation/NCBI_submission/zizania_palustris_13Nov2018_okGsv_renamedNCBI2.fasta'

for i in $(cat downsampling_2fold.txt); do
STEM=$(echo ${i} | cut -f 1 -d "/")
bwa mem -t 32 $FASTA ${i} 2> $STEM/${STEM}_downsampling_2fold_bwa.err | samtools sort -o $STEM/${STEM}_downsampling_2fold_sorted.bam 2> $STEM/${STEM}_downsampling_2fold_samtools_sort.err;
done
