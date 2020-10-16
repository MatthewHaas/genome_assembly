#!/bin/bash -l
#PBS -l nodes=1:ppn=8,mem=30g,walltime=6:00:00
#PBS -m abe
#PBS -M haasx092@umn.edu
#PBS -e run_jcvi.err
#PBS -o run_jcvi.out
#PBS -N run_jcvi

cd /home/jkimball/haasx092/other_synteny_figures/synteny_figure

export PATH=$PATH:/home/jkimball/haasx092/synteny_figure/last-1060/src
export PATH=$PATH:/home/jkimball/haasx092/synteny_figure/install-tl-20200505/1/bin/x86_64-linux

module load python

python -m jcvi.formats.gff bed rice.gene_structures_post_PASA_updates.21917.gff3 -o wild_rice.bed
python -m jcvi.formats.gff bed Osativa_323_v7.0.gene.gff3.gz -o oryza.bed
#python -m jcvi.formats.gff bed Oryza_rufipogon.OR_W1943.38.gff3.gz  -o  rufipogon.bed
python -m jcvi.formats.gff bed Zlat_genome_V1.gff.gz -o latifolia.bed # Careful with this BED file; it caused trouble in this pipeline

python -m jcvi.formats.fasta format rice.gene_structures_post_PASA_updates.nucleotides.fsa wild_rice.cds
python -m jcvi.formats.fasta format Osativa_323_v7.0.cds.fa.gz oryza.cds
python -m jcvi.formats.fasta format Oryza_rufipogon.OR_W1943.cds.all.fa.gz rufipogon.cds
python -m jcvi.formats.fasta format Zlat_V1.cds.fa.gz latifolia.cds

# The purpose of this section was to get the BED and CDS files to agree with one other in terms of chromosome and gene names
# Oryza rufipogon is included here because we attempted to include it in the figure, but synteny with O. sativa is almost perfect so we ultimately decided not to include it.
sed -i  's/\.MSUv7.0//g' oryza.bed
sed -i  's/\..*$//g' oryza.cds
sed -i 's/\-.*$//g' wild_rice.cds
#sed -i 's/\..*$//g' rufipogon.cds
#sed -i 's/gene://g' rufipogon.bed
sed -i 's/transcript_id//g' latifolia.bed
sed -i 's/"//g' latifolia.bed
sed -i 's/locus=.*//g' latifolia.cds
sed -i 's/ID=//g' latifolia.bed
sed -i 's/Parent=//g' latifolia.bed
sed -i 's/scaffold_/Chr/g' latifolia.bed

python -m jcvi.compara.catalog ortholog wild_rice oryza --cscore=1.0 --no_strip_names
python -m jcvi.compara.catalog ortholog oryza rufipogon --cscore=0.99 --no_strip_names
python -m jcvi.compara.catalog ortholog latifolia wild_rice --cscore=0.99 --no_strip_names

python -m jcvi.graphics.dotplot wild_rice.oryza.anchors
python -m jcvi.graphics.dotplot oryza.rufipogon.anchors
python -m jcvi.graphics.dotplot latifolia.wild_rice.anchors

python -m jcvi.compara.synteny screen --minspan=30 --simple wild_rice.oryza.anchors wild_rice.oryza.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple oryza.rufipogon.anchors oryza.rufipogon.anchors.new
python -m jcvi.compara.synteny screen --minspan=30 --simple latifolia.wild_rice.anchors latifolia.wild_rice.anchors.new

python -m jcvi.graphics.karyotype seqids layout --format png