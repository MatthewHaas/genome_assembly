#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=60g,walltime=6:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e clustalw.err
#PBS -o clustalw.out
#PBS -N clustalw


cd /home/jkimball/haasx092

module load clustalw

clustalw < clustalw_input.txt
