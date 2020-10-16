#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=15g,walltime=1:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e replace_scaffold_names.err
#PBS -o replace_scaffold_names.out
#PBS -N replace_scaffold_names
 
# Purpose of this code is to search specific repeat element density files for actual scaffold-specific numbers (from Dovetail release)
# and replace them with 1-16 numbering format in order to work with circos settings.
cd /home/jkimball/haasx092/circos/200511_version

for repeat_class in Gypsy Copia otherLTRs;
do
sed -i "s/zp1/zp1/g" ${repeat_class}/Scaffold_1_${repeat_class}_binned.txt # Not really necessary--replacing zp1 with zp1, but keeping for consistency with other scaffolds.
sed -i "s/zp3/zp2/g" ${repeat_class}/Scaffold_3_${repeat_class}_binned.txt
sed -i "s/zp7/zp3/g" ${repeat_class}/Scaffold_7_${repeat_class}_binned.txt
sed -i "s/zp9/zp4/g" ${repeat_class}/Scaffold_9_${repeat_class}_binned.txt
sed -i "s/zp13/zp5/g" ${repeat_class}/Scaffold_13_${repeat_class}_binned.txt
sed -i "s/zp18/zp6/g" ${repeat_class}/Scaffold_18_${repeat_class}_binned.txt
sed -i "s/zp48/zp7/g" ${repeat_class}/Scaffold_48_${repeat_class}_binned.txt
sed -i "s/zp51/zp8/g" ${repeat_class}/Scaffold_51_${repeat_class}_binned.txt
sed -i "s/zp70/zp9/g" ${repeat_class}/Scaffold_70_${repeat_class}_binned.txt
sed -i "s/zp93/zp10/g" ${repeat_class}/Scaffold_93_${repeat_class}_binned.txt
sed -i "s/zp415/zp11/g" ${repeat_class}/Scaffold_415_${repeat_class}_binned.txt
sed -i "s/zp693/zp12/g" ${repeat_class}/Scaffold_693_${repeat_class}_binned.txt
sed -i "s/zp1062/zp13/g" ${repeat_class}/Scaffold_1062_${repeat_class}_binned.txt
sed -i "s/zp1063/zp14/g" ${repeat_class}/Scaffold_1063_${repeat_class}_binned.txt
sed -i "s/zp1064/zp15/g" ${repeat_class}/Scaffold_1064_${repeat_class}_binned.txt
sed -i "s/zp1065/zp16/g" ${repeat_class}/Scaffold_1065_${repeat_class}_binned.txt
done