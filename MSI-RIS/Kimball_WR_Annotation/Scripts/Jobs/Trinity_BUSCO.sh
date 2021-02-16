#/bin/bash -l
#PBS -l nodes=1:ppn=12,mem=31gb,walltime=12:00:00
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -A jkimball
#PBS -W group_list=jkimball

module load python3/3.6.3_anaconda5.0.1
. /panfs/roc/msisoft/anaconda/anaconda3-5.0.1/etc/profile.d/conda.sh
conda activate busco4_msi

TRINITY_ASM="/home/jkimball/shared/WR_Annotation/Transcriptome_Asm_trinity/Trinity.fasta"
OUT="/home/jkimball/shared/WR_Annotation/Transcriptome_BUSCO"

cd "${OUT}"

busco \
    -m transcriptome \
    -i "${TRINITY_ASM}" \
    -o busco \
    -l poales_odb10 \
    -c 12
