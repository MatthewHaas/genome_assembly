
#!/bin/bash
#PBS -l mem=251gb,nodes=1:ppn=8,walltime=24:00:00
#PBS -A jkimball
#PBS -W group_list=jkimball
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q ram256g

# Load modules
module load hisat2/2.1.0
module load samtools/1.9

cd /home/jkimball/shared/WR_Annotation/HISAT2_Idx
hisat2-build -p 8 --ss ./splice_sites.txt genome.fa WR_genome
