#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account jkimball
#SBATCH -o get_pacbio_seqids.out
#SBATCH -e get_pacbio_seqids.err

cd /home/jkimball/haasx092

FILE_DIR='/home/jkimball/shared/PacBio_FilesJenny'

zcat ${FILE_DIR}/DTG-DNA-358_cell1.fastq.gz | grep "@" > cell1_headers.txt
zcat ${FILE_DIR}/DTG-DNA-358_cell2.fastq.gz | grep "@" > cell2_headers.txt
zcat ${FILE_DIR}/DTG-DNA-358_cell3.fastq.gz | grep "@" > cell3_headers.txt
zcat ${FILE_DIR}/DTG-DNA-358_cell4.fastq.gz | grep "@" > cell4_headers.txt
zcat ${FILE_DIR}/DTG-DNA-358_cell5.fastq.gz | grep "@" > cell5_headers.txt
zcat ${FILE_DIR}/DTG-DNA-358_cell6.fastq.gz | grep "@" > cell6_headers.txt
zcat ${FILE_DIR}/DTG-DNA-358_cell7.fastq.gz | grep "@" > cell7_headers.txt
zcat ${FILE_DIR}/DTG-DNA-358_cell8.fastq.gz | grep "@" > cell8_headers.txt

cat cell1_headers.txt cell2_headers.txt cell3_headers.txt cell4_headers.txt cell5_headers.txt cell6_headers.txt cell7_headers.txt cell8_headers.txt > all_headers_concat.txt
