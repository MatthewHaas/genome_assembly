-bash: conda: command not found
#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e filter_phylip_seqfiles.err
#PBS -o filter_phylip_seqfiles.out
#PBS -N filter_phylip_seqfiles

cd /home/jkimball/haasx092/duplicated_orthogroups

module load python3

for orthogroup in $(cat first_orthogroup_block.txt);
do
python filter_phylip_seqfiles.py $orthogroup/${orthogroup}_slim_backtranslated_with_space.phy $orthogroup/${orthogroup}_one_line_per_species_back.phy
done