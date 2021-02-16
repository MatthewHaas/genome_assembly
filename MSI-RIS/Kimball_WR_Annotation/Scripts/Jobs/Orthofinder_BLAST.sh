#!/bin/bash
#PBS -l pmem=1950mb,nodes=1:ppn=128,walltime=72:00:00
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q mangi
#PBS -A jkimball
#PBS -W group_list=jkimball

module load parallel
module load ncbi_blast+/2.8.1

# Run the BLAST searches separately
IN_CMDS="/home/jkimball/konox006/Projects/Kimball_WR_Annotation/Scripts/Jobs/Orthofinder_BLAST_Cmds_2020-06.txt"
# Run the BLAST searches in batches of 16
parallel -j 16 < "${IN_CMDS}"
