# 9 June 2020
# WD: /home/jkimball/haasx092/circos/200511_Version
# This R code is for getting coverage information from RNA-seq data to add to circos plot

# Load required packages
library(data.table)

# Read in data
x <- fread("/home/jkimball/shared/WR_Annotation/Tissue_Specificity/CHURP_Out/Counts/subread_counts.txt", header=TRUE)

# Simplify the Chr, Start, End, and Strand columns so that there is only one value for each
# Convert Start and End columns to numeric.
x[, Chr := sub(";.*$", "", Chr)]
x[, Start := as.numeric(sub(";.*$", "", Start))]
x[, End := as.numeric(sub(";.*$", "", End))]
x[, Strand := sub(";.*$", "", Strand)]

# Do the same for tissue-specific counts and convert to numeric
x[, Female := as.numeric(sub(";.*$", "", Female))]
x[, Flower := as.numeric(sub(";.*$", "", Flower))]
x[, Leaf := as.numeric(sub(";.*$", "", Leaf))]
x[, Male := as.numeric(sub(";.*$", "", Male))]
x[, Root := as.numeric(sub(";.*$", "", Root))]
x[, Sheath := as.numeric(sub(";.*$", "", Sheath))]
x[, Stem := as.numeric(sub(";.*$", "", Stem))]

# Extract gene ids from table
genes <- x$Geneid

# Define scaffolds of interest
scaffolds <- c("Scaffold_1", "Scaffold_3", "Scaffold_7", "Scaffold_9", "Scaffold_13", "Scaffold_18", 
			   "Scaffold_48", "Scaffold_51", "Scaffold_70", "Scaffold_93", "Scaffold_415",
			   "Scaffold_693", "Scaffold_1062", "Scaffold_1063", "Scaffold_1064", "Scaffold_1065")

# Keep only scaffolds of interest
x <- x[Chr %in% scaffolds]

# Find the mean read count for all tissues
# This will take a long time--there is probably a more efficient way at doing this.
for (i in genes) {
x[Geneid == i, average := mean(x[Geneid == i,c(Female, Flower, Leaf, Male, Root, Seed, Sheath, Stem)])]
}

# Keep only necessary columns
rna_coverage <- x[, .(Geneid, Chr, Start, average)]


# Find SNP density (almost identical to how retrotransposon density was found)
lapply(c(1,3,7,9,13,18,48,51,70,93,415,693,1062,1063,1064,1065), function (i) {
rna_coverage[Chr == paste0("Scaffold_", i), .(.N, average), key=.(bin=Start %/% 1000000 * 1000000)] -> rna

rna[, start := bin]
rna[, end := start + 1000000]
options(scipen = 999) # Gets rid of scientific notation

setcolorder(rna, c("start", "end", "N", "bin"))
rna[, bin := NULL]
}) -> rna_binned

# Average the RNA counts by 1 megabase bins
for(i in c(1:16)){
	for(j in seq(from = 0, to=max(rna_binned[[i]]$start), by=1e6)){
rna_binned[[i]][start == j, averaged_by_bin := mean(average)]
	}
}

# Keep unique values only
for (i in c(1:16)){
rna_binned[[i]] <- unique(rna_binned[[i]][, .(start, end, N, averaged_by_bin)])
}

# Add chromosome column, assign appropriate "zp" number, and remove the "N" column--no longer needed here
for (i in c(1:16)) {
rna_binned[[i]][, chr := paste0("zp", i)]
rna_binned[[i]][, N := NULL]
setcolorder(rna_binned[[i]], c("chr", "start", "end", "averaged_by_bin"))
}

# Name each item in the list
names(rna_binned) <- c("scf1", "scf3", "scf7", "scf9", "scf13", "scf18", "scf48", "scf51", "scf70", "scf93", "scf415", "scf693", "scf1062", "scf1063", "scf1064", "scf1065")

# Write each item to a separate text file
for (i in c(1:16)) {
write.table(rna_binned[[i]], file = paste0("scf", i, "_rna_counts"), sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

# Save data
save(x, rna_binned, file = "200609_rna_coverage.Rata")