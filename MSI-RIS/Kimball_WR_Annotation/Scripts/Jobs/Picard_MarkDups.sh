#!/bin/bash -l
#PBS -l nodes=1:ppn=1,mem=16gb,walltime=24:00:00
#PBS -M konox006@umn.edu
#PBS -m abe
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -q mangi

module load picard/2.18.16

#cd /home/jkimball/shared/WR_Annotation/Tissue_Specificity/Alignments
cd /panfs/roc/scratch/konox006/JKimball/RNASeq/Tissue_Specificity

BAM_LIST=($(find . -mindepth 1 -maxdepth 1 -name '*CoordSort.bam' -exec basename {} \; | sort -V))
CURR_BAM=${BAM_LIST[${PBS_ARRAYID}]}
_JAVA_OPTIONS="-Xmx15g" ${PTOOL}/picard.jar MarkDuplicates \
    I="${CURR_BAM}" \
    O="${CURR_BAM/.bam/_MarkDup.bam}" \
    REMOVE_DUPLICATES="false" \
    ASSUME_SORT_ORDER="coordinate" \
    M=$(basename "${CURR_BAM}" | cut -f 1 -d '_' )_Metrics.txt
