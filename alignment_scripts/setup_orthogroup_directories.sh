#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e setup_orthogroup_dirs.err
#PBS -o setup_orthogroup_dirs.out
#PBS -N setup_orthogroup_dirs

cd /home/jkimball/haasx092/duplicated_orthogroups

for orthogroup in $(cat first_orthogroup_block.txt);
do
mkdir ${orthogroup}
cp /home/jkimball/shared/WR_Annotation/Orthofinder_2020-06/Grass_Proteomes/OrthoFinder/Results_Jun02/WorkingDirectory/OrthoFinder/Results_Jun18/Orthogroup_Sequences/${orthogroup}.fa ./${orthogroup}
done
