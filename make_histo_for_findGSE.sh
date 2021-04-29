#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem=60g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o make_histo_for_findGSE.out
#SBATCH -e make_histo_for_findGSE.err

cd /home/jkimball/haasx092/findGSE

zcat *.fastq.gz | jellyfish count /dev/fd/0 -C -o test_21mer -m 21 -t 1 -s 5G
jellyfish histo -h 3000000 -o test_21mer.histo test_21mer
