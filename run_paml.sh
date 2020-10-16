#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=22g,walltime=24:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e paml.err
#PBS -o paml.out
#PBS -N paml

cd /home/jkimball/haasx092/Ks_substitution_rates

module load paml

yn00 runmode= -2, CodonFreq = 2 in codeml.ctl

## The following lines would go into the configuration (ctl) file;
## These are not the names of my files; I this from the documentation and am trying to figure out what they do
seqfile = abglobin.nuc # sequence data file name
outfile = substitution_test_out # main result file
verbose = 0 # 1: detailed output (list sequences), 0: concise output
icode = 0 # 0:universal code; 1:mammalian mt; 2-10:see below weighting = 0 * weighting pathways between codons (0/1)?
commonf3x4 = 0