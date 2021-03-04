#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=4:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=haasx092@umn.edu
#SBATCH -p amdsmall
#SBATCH --account=jkimball
#SBATCH -o scythe_mpileup_4fold.out
#SBATCH -e scythe_mpileup_4fold.err

cd /home/jkimball/haasx092/pilot_GBS

module load samtools
module load bcftools
module load htslib
module load parallel

export bams="downsampled_4fold_bam_list"
export prefix="200519_downsampled_4fold"
export ref="/home/jkimball/shared/WR_Annotation/NCBI_submission/zizania_palustris_13Nov2018_okGsv_renamedNCBI2.fasta"

parallel_samtools_processes=15

mkdir -p $prefix
cut -f 1 ${ref}.fai


scythe_mpileup() {
	REGIONS=${1}
	SCAFFOLD=$(echo ${REGIONS} | cut -f 1 -d ";")
	samtools mpileup -q 20 -gDVu \
		-b $bams \
		-r ${REGIONS} \
		-f ${ref} \
		| bcftools call -mv \
		| bgzip -c \
		> $prefix/${prefix}_${SCAFFOLD}.vcf.gz \
		2> $prefix/${prefix}_${SCAFFOLD}.err
}
  
export -f scythe_mpileup
  
parallel --will-cite -j $parallel_samtools_processes scythe_mpileup ::: $(cut -f 1 ${ref}.fai)
