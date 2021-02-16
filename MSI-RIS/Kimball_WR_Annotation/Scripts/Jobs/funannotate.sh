#!/bin/bash
#PBS -l mem=62GB,nodes=1:ppn=24,walltime=250:00:00
#PBS -A jkimball
#PBS -m abe
#PBS -j oe
#PBS -M mmacchie@umn.edu
#PBS -q max
#PBS -W group_list=jkimball
#PBS -N funannotate_update

module load funannotate/1.5.1
#module unload gmap
#module load gmap/2012-10-31

export EVM_HOME=/panfs/roc/msisoft/evidencemodeler/1.1.1/
export AUGUSTUS_CONFIG_PATH=/scratch.global/mmacchie/config/
export AUGUSTUS_BIN_PATH=/scratch.global/mmacchie/bin/
export GENEMARK_PATH=/panfs/roc/msisoft/genemark/4.32/
export BAMTOOLS_PATH=/panfs/roc/msisoft/bamtools/2.5.1/bin/
export PATH=/panfs/roc/msisoft/augustus/3.2.3.CentOS7/scripts/:$PATH

echo -n "Ran: "
date

cd /home/jkimball/shared/WR_Annotation/Annotation/ 
GENOME='/home/jkimball/shared/WR_Annotation/genome/reformatted_genome/zizania_palustris_13Nov2018_okGsv.fasta.masked.headerreformatted'
TRINITY='/home/jkimball/shared/WR_Annotation/Transcriptome_Asm_trinity/Trinity.fasta'
RNABAM='/home/jkimball/shared/WR_Annotation/Aligned_reads/rice_rnas_chromreformattedAligned.sortedByCoord.out.bam'
FASTQ_R1='/home/jkimball/shared/WR_Annotation/Trimmed_reads/WR_8Tissues.R1.trimmed.fastq'
FASTQ_R2='/home/jkimball/shared/WR_Annotation/Trimmed_reads/WR_8Tissues.R2.trimmed.fastq'

#--- check species options
# funannotate species

#---Predict gene models ----#
# run on amdsmall; mem=62GB,nodes=1:ppn=32,walltime=96:00:00
#funannotate predict -i ${GENOME} --species "rice" --transcript_evidence ${TRINITY} --rna_bam ${RNABAM} -o wild_rice3 --cpus 32 
# "--species or --augustus_species?"

# run on max; mem=62GB,nodes=1:ppn=24,walltime=120:00:00
funannotate update -i wild_rice3 --cpus 24 --left ${FASTQ_R1} --right ${FASTQ_R2} --trinity ${TRINITY} --jaccard_clip

# run on max; mem=2GB,nodes=1:ppn=1,walltime=300:00:00
# funnannotate update threw an error at the PASA step - required 1 core for the sqlite database. So I am going to run the step that crashed with one core, and try to resume the job.
#cd /home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/pasa/
#/panfs/roc/msisoft/pasa/2.3.3/scripts/cDNA_annotation_comparer.dbi -G /panfs/roc/groups/1/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/genome.fa --CPU 1 -M '/panfs/roc/groups/1/jkimball/shared/WR_Annotation/Annotation/rice'  > pasa_run.log.dir/rice.annotation_compare.32323.out

# run on 2020-01-24
#mkdir /home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/pasa2/
#cd /home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/pasa2/
#/panfs/roc/msisoft/pasa/2.3.3/scripts/cDNA_annotation_comparer.dbi -G /panfs/roc/groups/1/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/genome.fa --CPU 24 -M '/panfs/roc/groups/1/jkimball/shared/WR_Annotation/Annotation/rice'  > pasa_run.log.dir/rice.annotation_compare.32323.out


echo -n "Done: "
date



