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
#SBATCH -o run_ltr_harvest.out
#SBATCH -e run_ltr_harvest.err

module load genometools

cd /home/jkimball/haasx092/LTR_assembly_index

gt ltrharvest -index zizania_palustris_13Nov2018_okGsv_renamedNCBI2_copy \
-similar 90 -vic 10 -seed 20 -seqids yes \ -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 \
-motif TGCA -motifmis 1 > genome.harvest.scn
