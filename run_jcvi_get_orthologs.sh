#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account jkimball
#SBATCH -o jcvi_get_orthologs.out
#SBATCH -e jcvi_get_orthologs.err

cd /home/jkimball/haasx092/ortholog_ratios

export PATH=$PATH:/home/jkimball/haasx092/synteny_figure/last-1060/src
export PATH=$PATH:/home/jkimball/haasx092/synteny_figure/install-tl-20200505/1/bin/x86_64-linux

module load python

python -m jcvi.formats.gff bed ~/synteny_figure/Osativa_323_v7.0.gene.gff3.gz -o oryza.bed
python -m jcvi.formats.gff bed Zlat_genome_V1.gff -o latifolia.bed


python -m jcvi.formats.fasta format ~/synteny_figure/Osativa_323_v7.0.cds.fa.gz oryza.cds
python -m jcvi.formats.fasta format Zlat_V1.cds.fa.gz latifolia.cds

sed -i  's/\.MSUv7.0//g' oryza.bed
sed -i  's/\..*$//g' oryza.cds

python -m jcvi.compara.catalog ortholog latifolia oryza --cscore=0.99 --no_strip_names

python -m jcvi.compara.synteny depth --histogram latifolia.oryza.anchors
