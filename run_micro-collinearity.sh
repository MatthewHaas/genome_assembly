#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=30g,walltime=6:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e micro-collinearity.err
#PBS -o micro-collinearity.out
#PBS -N micro-collinearity

cd /home/jkimball/haasx092/other_synteny_figures/synteny_figure

export PATH=$PATH:/home/jkimball/haasx092/synteny_figure/last-1060/src
export PATH=$PATH:/home/jkimball/haasx092/synteny_figure/install-tl-20200505/1/bin/x86_64-linux

module load python

cat wild_rice.bed oryza.bed > wild_rice_oryza.bed

python -m jcvi.compara.synteny mcscan wild_rice.bed wild_rice.oryza.lifted.anchors --iter=1 -o wild_rice.oryza.i1.blocks

python -m jcvi.graphics.synteny blocks wild_rice_oryza.bed blocks.layout
