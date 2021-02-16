#!/bin/bash
#PBS -l mem=16gb,nodes=1:ppn=6,walltime=4:00:00
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q mangi

# Load modules
module load hisat2/2.1.0
module load samtools/1.9

# Set directories
READS_DIR="/panfs/roc/scratch/konox006/JKimball/RNASeq/Trimmed_Reads"
OUT_DIR="/home/jkimball/shared/WR_Annotation/Tissue_Specificity/Alignments"
GENOME="/home/jkimball/shared/WR_Annotation/HISAT2_Idx/WR_genome"

# Find all R1 files in the reads directory. Use their names to build the R2
R1=($(find "${READS_DIR}" -mindepth 1 -maxdepth 1 -type f -name '*1P.fq.gz' | sort -V))
CURR_R1="${R1[${PBS_ARRAYID}]}"
R2="${CURR_R1/1P/2P}"

# Build the samplename from the read name
SAMPLENM=$(basename "${CURR_R1}" | cut -f 1 -d '_')

# Align the reads with HISAT2
hisat2 \
    -p "${PBS_NUM_PPN}" \
    --no-mixed \
    --new-summary \
    --rna-strandness RF \
    -x "${GENOME}" \
    -1 "${CURR_R1}" \
    -2 "${R2}" \
    2> "${OUT_DIR}/${SAMPLENM}_Aln_Summary.txt" \
    | samtools view -hb -o "${OUT_DIR}/${SAMPLENM}.bam" -
