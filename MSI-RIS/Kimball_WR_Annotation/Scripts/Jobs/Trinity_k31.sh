#!/bin/bash
#PBS -l mem=998gb,nodes=1:ppn=32,walltime=24:00:00
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q ram1t

# Load modules
module load java/jdk1.8.0_171
module load bowtie2/2.3.4.1.CentOS7
module load samtools/1.7
module load R/3.5.0
module load ncbi_blast+/2.7.1.CentOS7
module load jellyfish/2.1.3
module load python3/3.6.3_anaconda5.0.1
module load salmon/0.14.1

# Define paths
TRINITY="/home/jkimball/shared/Software/trinityrnaseq-v2.8.6/Trinity"
READS_DIR="/panfs/roc/scratch/konox006/JKimball/RNASeq/Trimmed_reads"
OUTPUT_DIR="/panfs/roc/scratch/konox006/JKimball/RNASeq/trinity_k31_Transcriptome_Asm"

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

# Set some more parameters for the assembly
MIN_CONTIG=250

# Run Trinity
${TRINITY} \
    --seqType fq \
    --max_memory 996G \
    --left "${READS_DIR}/WR_8Tissues.R1.trimmed.fastq.gz" \
    --right "${READS_DIR}/WR_8Tissues.R2.trimmed.fastq.gz" \
    --SS_lib_type RF \
    --CPU 32 \
    --KMER_SIZE 31\
    --monitoring \
    --monitor_sec 30 \
    --min_contig_length ${MIN_CONTIG} \
    --output "${OUTPUT_DIR}"
