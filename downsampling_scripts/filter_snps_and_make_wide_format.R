# Load required libraries
library(data.table)
library(reshape2)
args <- commandArgs(trailingOnly = TRUE)
# Read in the data, give column names, and remove a unnecessary/empty column
fread(args[1]) -> x
setnames(x, c("scaffold", "position", "ref", "alt", "quality", "sample", "GT", "V8", "DP", "DV"))
x[, V8 := NULL]
# Clean up sample names by stripping the "relative path" style, leaving only the sample name.
# The sample names will be turned into column names for part of the analysis, so this will make it cleaner.
x[, sample := sub("/.+$", "", sample)]
x[, scaffold := sub(";.+$", "", scaffold)]
# Name scaffolds of interest
scaffolds_of_interest = c("ZPchr0001", "ZPchr0002", "ZPchr0003", "ZPchr0004", "ZPchr0005", "ZPchr0006", "ZPchr0007", "ZPchr0008", "ZPchr0009",
			  "ZPchr0010", "ZPchr0011", "ZPchr0012", "ZPchr0013", "ZPchr0014", "ZPchr0015", "ZPchr0016", "ZPchr0458")
# Retain only scaffolds of interest
x[scaffold %in% scaffolds_of_interest] -> y
# Filter based on genotype. Keep all homozygous calls but filter out heterozygotes with less than 3 ALT alleles
y -> z # removed filtering step, so adding this so the object names still work
# Convert from long format to wide format
dcast(z, scaffold + position + ref + alt ~ sample, value.var="GT") -> zz
zz <- as.data.table(zz)
# Write the filtered SNP table to a file
write.csv(zz, file=args[2], row.names=FALSE, col.names=TRUE)
# Save data (especially the most important ones--so I can retrieve lower depth--5 to 9 reads-- if necessary)
save(zz, y, scaffolds_of_interest, file=args[3])
