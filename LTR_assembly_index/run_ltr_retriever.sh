#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=10
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball 
#SBATCH -o run_ltr_retriever.out
#SBATCH -e run_ltr_retriever.err

module load ncbi_blast+/2.8.1
module load cd-hit
module load hmmer
module load repeatmasker

cd /home/jkimball/haasx092/LTR_assembly_index

LTR_retriever/LTR_retriever -genome zizania_palustris_13Nov2018_okGsv_renamedNCBI2_copy.fasta -inharvest genome.harvest.scn
