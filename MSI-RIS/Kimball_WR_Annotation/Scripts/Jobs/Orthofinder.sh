#!/bin/bash
#PBS -l mem=998gb,nodes=1:ppn=32,walltime=96:00:00
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q ram1t
#PBS -A jkimball
#PBS -W group_list=jkimball

module load python3/3.6.3_anaconda5.0.1
module load ncbi_blast+/2.8.1

# Set paths
#ORTHOFINDER="/home/jkimball/shared/Software/OrthoFinder_2.3.11/OrthoFinder/orthofinder"
ORTHOFINDER="/home/jkimball/shared/Software/OrthoFinder_Source_2.3.11/OrthoFinder_source/orthofinder.py"
BLAST_DIR="/home/jkimball/shared/WR_Annotation/Orthofinder_2020-06/Grass_Proteomes/OrthoFinder/Results_Jun02/WorkingDirectory"

python "${ORTHOFINDER}" \
    -a "${PBS_NUM_PPN}" \
    -b "${BLAST_DIR}"
