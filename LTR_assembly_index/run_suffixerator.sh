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
#SBATCH -o run_suffixerator.out
#SBATCH -e run_suffixerator.err

module load genometools

cd /home/jkimball/haasx092/LTR_assembly_index

gt suffixerator -db zizania_palustris_13Nov2018_okGsv_renamedNCBI2_copy.fasta -indexname zizania_palustris_13Nov2018_okGsv_renamedNCBI2_copy -tis -suf -lcp -des -ssp -sds -dna 2> suffixerator.err
