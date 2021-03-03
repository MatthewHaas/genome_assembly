# Load required packages
library(data.table)

# Set directory
setwd("~/Documents/wild_rice")

# Read in data
x <- fread("210302_samtools_ZPchrs_filtered_concat.csv")

# Make data.table
shattering_genes <- data.table(gene = c("ZPchr0010g10516", "ZPchr0010g7757",
                                      "ZPchr0003g18426", "ZPchr0013g34051",
                                      "ZPchr0458g22499", "ZPchr0001g31104",
                                      "ZPchr0005g15825", "ZPchr0010g10516",
                                      "ZPchr0010g7757", "ZPchr0002g26578", 
                                      "ZPchr0004g39486", "ZPchr0006g44379",
                                      "ZPchr0006g42764", "ZPchr0228g22246",
                                      "ZPchr0006g44581"),
                               chrom = c("ZPchr0010", "ZPchr0010",
                                         "ZPchr0003", "ZPchr0013",
                                         "ZPchr0458", "ZPchr0001",
                                         "ZPchr0005", "ZPchr0010",
                                         "ZPchr0010", "ZPchr0002",
                                         "ZPchr0004", "ZPchr0006",
                                         "ZPchr0006", "ZPchr0228",
                                         "ZPchr0006"),
                               start = c(50109711, 50310476, 53532217,
                                         67454025, 3346978, 13965773,
                                         10456511, 50109711, 50310476,
                                         29878612, 97288584, 53321051,
                                         53486138, 8170, 28201679))

# Count SNPs within 1 MB in both directions of the shattering genes listed ni Table 3
for (i in shattering_genes$gene){
  shattering_genes[gene == i, snp_count := nrow(x[scaffold == shattering_genes[gene == i]$chrom & 
                                                                    (position >= shattering_genes[gene == i]$start - 1e6) & 
                                                                    (position <= shattering_genes[gene == i]$start + 1e6)])]
}    
# Write to a file
write.csv(shattering_genes, file = "210302_snp_count_for_shattering_genes.csv")
