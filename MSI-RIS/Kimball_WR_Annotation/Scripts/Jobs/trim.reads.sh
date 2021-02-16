#!/bin/bash
#PBS -l mem=4GB,nodes=1:ppn=2,walltime=12:00:00
#PBS -A jkimball 
#PBS -m abe
#PBS -j oe
#PBS -M mmacchie@umn.edu
#PBS -q small
#PBS -W group_list=jkimball
#PBS -N trimmomatic

module load trimmomatic/0.33
module load parallel/20160822

READS_DIR='/panfs/roc/scratch/konox006/JKimball/RNASeq/Cleaned_Reads'


echo -n "Ran: "
date

# define the trimmomatic function
trimmo() {
	FQ1=${1}
	SAMP=${FQ1/_R1.fastq.gz/}
	FQ2=${SAMP}_R2.fastq.gz
	READS_DIR2='/panfs/roc/scratch/konox006/JKimball/RNASeq/Trimmed_reads'
	java -jar $TRIMMOMATIC/trimmomatic.jar PE -threads 2 ${FQ1} ${FQ2} ${READS_DIR2}/${SAMP}.R1.trimmed.fastq.gz ${READS_DIR2}/${SAMP}.R1.SE.trimmed.fastq.gz ${READS_DIR2}/${SAMP}.R2.trimmed.fastq.gz ${
READS_DIR2}/${SAMP}.R2.SE.trimmed.fastq.gz ILLUMINACLIP:$TRIMMOMATIC/adapters/all_illumina_adapters.fa:2:30:10:2:true 
}

export -f trimmo

cd ${READS_DIR}

parallel --jobs 2 --joblog parallel_trimming_logfile.txt trimmo ::: *_R1.fastq.gz

echo -n "Done: "
date

