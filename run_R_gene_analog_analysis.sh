#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=50gb,walltime=48:00:00
#PBS -q amdsmall
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -N RGA_analysis
#PBS -e RGA_analysis.err
#PBS -o RGA_analysis.out

# This script is for conducting Resistance Gene Analog (RGA) analysis for Northern Wild Rice
# Uses the RGAugury pipeline (Li, P., Quan, X., Jia, G., Xiao, J., Cloutier, S. and You, F.M. (2016) RGAugury: a pipeline for genome-wide prediction of resistance gene analogs (RGAs) in plants. BMC genomics, 17, 852.)

# Change into proper directory
# NOTE: do not leave as variables, change to your group and username
cd /home/jkimball/haasx092/RGA

# Get correct version of gcc
module load gcc/5.1.0
export LD_PRELOAD=/panfs/roc/msisoft/gcc/5.1.0/lib64/libgomp.so
export LD_PRELOAD=/panfs/roc/msisoft/gcc/5.1.0/lib64/libstdc++.so:$LD_PRELOAD
export LD_PRELOAD=/panfs/roc/msisoft/gcc/5.1.0/lib64/libstdc++.so.6:$LD_PRELOAD

# Load required software
module load bioperl
module load hmmer
module load phobius/1.01
module load ncbi_blast+
module load interproscan/5.23-62.0

# Load libraries
module load zlib
module load freetype
module load libgd

# Set up environment
# NOTE: each reference to $GROUP and $USER should again be changed to your group and username
### LIBPNG
export PATH=$PATH:/home/jkimball/haasx092/RGA/bin/libpng-1.6.37/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jkimball/haasx092/RGA/bin/libpng-1.6.37/lib
### pfam_scan.pl
export PATH=$PATH:/home/jkimball/haasx092/RGA/bin/PfamScan
### RGAugury
export PATH=$PATH:/home/jkimball/haasx092/RGA/bin/rgaugury
### cCoils
export COILSDIR=/home/jkimball/haasx092/RGA/bin/rgaugury/coils # path to scoils-ht
export PATH=$PATH:$COILSDIR # path to scoils-ht
### CViT
export PATH=$PATH:/home/jkimball/haasx092/RGA/bin/cvit.1.2.1 
### PfamDB
export PFAMDB=/home/jkimball/haasx092/RGA/bin/pfamdb 
export PERL5LIB=/home/jkimball/haasx092/RGA/bin/PfamScan

# Define filenames
PROTEIN_FASTA=/home/jkimball/haasx092/RGA/inputFiles/augustus.pasa.proteins.fa
GENOME_FASTA=/home/jkimball/haasx092/inputFiles/Zizania_palustris_0.1_CDS.fasta.gz
GFF_FILE=/home/jkimball/haasx092/RGA/inputFiles/rice.gene_structures_post_PASA_updates.21917.gff3 

# Other
PREFIX=ZizaniaPalustris
# It is also possible to specify the cpu or threads number. For now, will leave as default (default=2).

# Run pipeline
perl bin/rgaugury/RGAugury.pl -p $PROTEIN_FASTA -g $GENOME_FASTA -gff $GFF_FILE -pfx $PREFIX
