#!/bin/bash -l
#PBS -l nodes=1:ppn=1,mem=2gb,walltime=36:00:00
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -M konox006@umn.edu
#PBS -m abe
#PBS -q mesabi

cd /home/jkimball/shared/WR_Annotation
tar -czvf WR_Trinity.tar.gz Transcriptome_Asm_trinity/
