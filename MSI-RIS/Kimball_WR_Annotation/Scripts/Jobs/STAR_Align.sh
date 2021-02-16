#!/bin/bash
#PBS -l mem=20gb,nodes=1:ppn=6,walltime=4:00:00
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q mangi

# Load modules
module load star/2.7.1a
module load samtools/1.9

# Set directories
READS_DIR="/panfs/roc/scratch/konox006/JKimball/RNASeq/Trimmed_Reads"
OUT_DIR="/panfs/roc/scratch/konox006/JKimball/RNASeq/Tissue_Specificity/STAR_Aln"
GENOME="/home/jkimball/shared/WR_Annotation/genome/reformatted_genome/"

# Find all R1 files in the reads directory. Use their names to build the R2
R1=($(find "${READS_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*1P.fq.gz' | sort -V))
CURR_R1="${R1[${PBS_ARRAYID}]}"
R2="${CURR_R1/1P/2P}"

# Build the samplename from the read name
SAMPLENM=$(basename "${CURR_R1}" | cut -f 1 -d '_')

STAR \
    --runMode alignReads \
    --runThreadN "${PBS_NUM_PPN}" \
    --genomeDir "${GENOME}" \
    --readFilesCommand 'zcat' \
    --readFilesIn "${CURR_R1}" "${R2}" \
    --outFileNamePrefix "${OUT_DIR}/${SAMPLENM}" \
    --limitBAMsortRAM 24846526061 \
    --outSAMtype BAM SortedByCoordinate
