#!/bin/bash -l
#PBS -l nodes=1:ppn=6,mem=16gb,walltime=12:00:00
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q mesabi

module load fastqc/0.11.8

IN_DIR="/home/jkimball/data_release/umgc/novaseq/190828_A00223_0198_BHCMFVDRXX/Kimball_Project_003"
OUT_DIR="/panfs/roc/scratch/konox006/JKimball/RNASeq/FastQC"

mkdir -p "${OUT_DIR}"
fastqc -t 6 -f fastq -o "${OUT_DIR}" "${IN_DIR}"/*.fastq.gz
