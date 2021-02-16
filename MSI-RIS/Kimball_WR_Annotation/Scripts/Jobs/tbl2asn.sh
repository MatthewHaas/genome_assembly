#!/bin/bash
#PBS -l mem=60GB,nodes=1:ppn=4,walltime=96:00:00
#PBS -A jkimball
#PBS -m abe
#PBS -j oe
#PBS -M mmacchie@umn.edu
#PBS -q amd2tb
#PBS -W group_list=jkimball
#PBS -N tbl2asn

cd /home/jkimball/shared/WR_Annotation/NCBI_submission 

#./linux64.tbl2asn -i ../genome/reformatted_genome/zizania_palustris_13Nov2018_okGsv.fasta.masked.headerreformatted.fa -f feature.table.test.tbl -t template.sbt -o zpal.asn -V vbt -s T -m B -l paired-ends -a r10k -W T

./linux64.tbl2asn -i ../genome/zizania_palustris_13Nov2018_okGsv_renamedNCBI2.fasta -f feature.table.NCBI.tbl -t ncbi_template.sbt -o zpal.asn -V vbt -s T -m B -l paired-ends -a r10k -W T -j "[organism=Zizania palustris]"






