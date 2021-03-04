# README for downsampling scripts

There are a lot of very similar files in this directory. Usually, I try to make general scripts so that they can be recycled for nearly identical scripts. The script [filter_snps_and_make_wide_format.R](filter_snps_and_make_wide_format.R) is a good example of this. In the case of downsampling, I used different subdirectories for some of the analyses so the only difference might be the path to the working directory. The input files might also differ slightly. For example, in [filter_with_vcftools_2fold.sh](filter_with_vcftools_2fold.sh) (as well as the 4-fold and 8-fold versions), the input file name [vcf_list_new_2fold.txt](vcf_list_new_2fold.txt) was just called ```vcf_list_new.txt``` in the directory, but has been renamed here so they won't conflict with each other.

# Workflow

1. Trimmed FASTQ files were downsampled using [downsampling_2fold.sh](downsampling_2fold.sh), [downsampling_4fold.sh](downsampling_4fold.sh), and [downsampling_8fold.sh](downsampling_8fold.sh).
2. Downsampled FASTQ files were aligned to the NWR genome using [run_bwa_downsampling_2fold.sh](run_bwa_downsampling_2fold.sh), [run_bwa_downsampling_4fold.sh](run_bwa_downsampling_4fold.sh), and [run_bwa_downsampling_8fold.sh](run_bwa_downsampling_8fold.sh).
3. BAM files were indexed using [make_bam_indicies.sh](make_bam_indicies.sh).
4. SNP calling was done (with samtools mpileup/bcftools call) using [scythe_mpileup_downsampled_2fold.sh](scythe_mpileup_downsampled_2fold.sh), [scythe_mpileup_downsampled_4fold.sh](scythe_mpileup_downsampled_4fold.sh), and [scythe_mpileup_downsampled_8fold.sh](scythe_mpileup_downsampled_8fold.sh).
5. VCF files were filtered (10% missing data, bi-allelic sites only, minimum of 6 reads) using [filter_with_vcftools_2fold.sh](filter_with_vcftools_2fold.sh), [filter_with_vcftools_4fold.sh](filter_with_vcftools_4fold.sh), and [filter_with_vcftools_8fold.sh](filter_with_vcftools_8fold.sh).
