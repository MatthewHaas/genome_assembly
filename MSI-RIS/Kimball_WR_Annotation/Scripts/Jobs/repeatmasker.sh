#!/bin/bash
#PBS -l mem=32GB,nodes=1:ppn=24,walltime=200:00:00
#PBS -A jkimball
#PBS -m abe
#PBS -j oe
#PBS -M mmacchie@umn.edu
#PBS -q max
#PBS -W group_list=jkimball
#PBS -N repeatmasker

module load repeatmodeler/1.0.11 repeatmasker/4.0.5 perl/modules.centos7.5.26.1


echo -n "Ran: "
date

GENOME='/home/jkimball/shared/WR_Annotation/genome/zizania_palustris_13Nov2018_okGsv.fasta'

#cd /panfs/roc/scratch/konox006/JKimball/RNASeq/genome/

# 1. Run RepeatModeler
## Build repeat database
#BuildDatabase -name oryza_sativa -engine ncbi ${GENOME}

#echo -n "Done: "
#date


## Given a genomic database from above, build, refine, and classify consensus models of interspersed repeats 
#RepeatModeler -engine ncbi -pa 24 -database oryza_sativa

## Combine repeat modeler results with repeatmasker (RM) database
### Convert RM.embl to fasta
#buildRMLibFromEMBL.pl /panfs/roc/msisoft/repeatmasker/4.0.5/Libraries/RepeatMaskerLib.embl > RepeatMaskerLib.fasta
### Combine RM fasta with RM fasta (consensi.fa.classified)
#cat RepeatMaskerLib.fasta consensi.fa.classified > combined_repeat_libs.fasta


cd /home/jkimball/shared/WR_Annotation/repeatmasker/
# 2. Run RepeatMasker
RepeatMasker -xsmall -x -gff -pa 6 -s -lib combined_repeat_libs.fasta -e ncbi ${GENOME}


echo -n "Done: "
date



