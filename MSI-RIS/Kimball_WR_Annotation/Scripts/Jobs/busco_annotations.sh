#/bin/bash -l
#PBS -l nodes=1:ppn=12,mem=31gb,walltime=12:00:00
#PBS -m abe
#PBS -M mmacchie@umn.edu
#PBS -A jkimball
#PBS -W group_list=jkimball


module load python3/3.6.3_anaconda5.0.1
. /panfs/roc/msisoft/anaconda/anaconda3-5.0.1/etc/profile.d/conda.sh
conda activate /home/msistaff/konox006/.conda/envs/busco4_msi/

PASA="/home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/update_misc/pasa/augustus.pasa.proteins.longestiso.fa"
OUT="/home/jkimball/shared/WR_Annotation/funannotate_BUSCO"

cd "${OUT}"

busco \
    -m proteins \
    -i "${PASA}" \
    -o busco \
    -l poales_odb10 \
    -c 12
