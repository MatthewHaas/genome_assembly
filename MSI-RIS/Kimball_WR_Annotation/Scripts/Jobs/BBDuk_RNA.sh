#!/bin/bash -l
#PBS -l nodes=1:ppn=12,mem=31gb,walltime=24:00:00
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q mesabi

module load java/openjdk-11.0.2

# Use the SILVA database to screen for rRNA contamination. This reference file
# was generated for CHURP!
SILVA="/home/msistaff/public/CHURP_Deps/v0/db/SILVA_132_LSU_SSU_Ref_Dedup_Kmers_min100.fasta.gz"
BBDUK="/home/jkimball/shared/Software/bbmap/bbduk.sh"

READS_DIR="/home/jkimball/data_release/umgc/novaseq/190828_A00223_0198_BHCMFVDRXX/Kimball_Project_003"
OUT_DIR="/panfs/roc/scratch/konox006/JKimball/RNASeq/rRNA_Screen"

# Make an array of the R1 sequences
R1=($(find "${READS_DIR}" -maxdepth 1 -mindepth 1 -type f -name '*R1_001.fastq.gz' | sort -V))

# Grab the current R1 using the array ID
CURR_R1="${R1[${PBS_ARRAYID}]}"
# Get the basename because we will need it to define output names
R1BNAME=$(basename ${CURR_R1})
R2BNAME="${R1BNAME/_R1_001.fastq.gz/_R2_001.fastq.gz}"

# run BBDUK
_JAVA_OPTIONS="-Xmx30g" "${BBDUK}" \
    in="${CURR_R1}" \
    in2="${READS_DIR}/${R2BNAME}" \
    ref="${SILVA}" \
    k=25 \
    editdistance=1 \
    out="${OUT_DIR}/${R1BNAME/_R1_001.fastq.gz/_R1_rRNADep.fastq.gz}" \
    out2="${OUT_DIR}/${R2BNAME/_R2_001.fastq.gz/_R2_rRNADep.fastq.gz}" \
    outm="${OUT_DIR}/${R1BNAME/_R1_001.fastq.gz/_R1_rRNAMatch.fastq.gz}" \
    outm2="${OUT_DIR}/${R2BNAME/_R2_001.fastq.gz/_R2_rRNAMatch.fastq.gz}" \
    outs="${OUT_DIR}/${R1BNAME/_R1_001.fastq.gz/_Single_rRNAMatch.fastq.gz}" \
    stats="${OUT_DIR}/${R1BNAME/_R1_001.fastq.gz/_BBDuk_Stats.txt}" \
    ordered=true \
    prealloc=t
