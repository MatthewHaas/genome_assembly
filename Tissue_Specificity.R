# This script written by Tom Kono
# Calculate a tissue specificity index (tau) for each gene annotated in the
# assembly. We might have to exclude the root sample due to ribosomal RNA
# contamination.

library(viridis)
library(pheatmap)
library(edgeR)

# Tissue specificity threshold
TAU_THRESHOLD <- 0.8

# Filtering criteria. We will skip the variance filter here because that is
# part of the signal
MIN_LEN <- 200
MIN_CTS <- 10

# Subread counts matrix
raw_mat <- read.table("/Users/tomkono/Dropbox/GitHub/RIS/Kimball_WR_Annotation/Results/RNASeq_Alignment/subread_counts.txt.gz", header=TRUE, sep="\t")

# Apply the filtering criteria to the matrix
# Length filtering is easy
keep_len <- raw_mat$Length >= MIN_LEN
# Define a function to do the counts filtering and apply it over the matrix. We
# are doing a twist on the counts filtering - at least X samples must have
# counts This is opposed to the DEG filtering criteria, where at most X samples
# must have 0 counts.
flt_cts <- function(counts, size, minexp=MIN_CTS) {
    pass <- as.numeric(counts) >= minexp
    if (sum(pass) <= size) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}
# Set the smallest group size to 1 because there is one replicate per tissue.
keep_cts <- apply(raw_mat[,seq(-1, -6)], 1, flt_cts, 1)
# those passing both criteria are kept
keep <- keep_len & keep_cts
flt_mat <- raw_mat[keep,]

# Print some summaries
sum(keep_len)
sum(keep_cts)
dim(raw_mat)
dim(flt_mat)

# Covert the filtered matrix into an EdgeR data object and calculate FPKM
edge_mat <- DGEList(counts=flt_mat[,seq(-1,-6)], genes=flt_mat[,1])
edge_mat <- calcNormFactors(edge_mat)
flt_rpkm <- rpkm(edge_mat, gene.length=flt_mat$Length)
flt_rpkm <- as.data.frame(flt_rpkm)

# Now we will calculate tau
calc_tau <- function(x) {
    expression <- as.numeric(x)
    xhat <- expression / max(expression)
    s <- 1-xhat
    tau <- sum(s)/(length(x) - 1)
    return(tau)
}
flt_rpkm$Tau <- as.numeric(apply(flt_rpkm, 1, calc_tau))
flt_rpkm$Gene <- as.character(raw_mat$Geneid[keep])

# Add a column for the genes that are specific to which tissue
get_tissue <- function(x) {
    expression <- as.numeric(x[1:(length(x)-2)])
    tissues <- names(x[1:(length(x)-2)])
    tau <- x["Tau"]
    if(tau < TAU_THRESHOLD) {
        return(NA)
    } else {
        return(tissues[which.max(expression)])
    }
}

# Write the tau values to disk
flt_rpkm$Tissue.Specific <- apply(flt_rpkm, 1, get_tissue)
write.csv(flt_rpkm, file="WR_Tau_AllTissues.csv", quote=FALSE, row.names=FALSE)

# Make a nice lil plot
pdf(file="Tau_Density_AllTissues.pdf", height=8, width=8)
tdens <- flt_rpkm$Tau[!is.na(flt_rpkm$Tau)]
hist(
    tdens,
    xlim=c(0, 1),
    breaks=100,
    col="lightgrey",
    border=NA,
    xlab="Tau",
    ylab="Count",
    main="Tissue Specificity in Northern Wild Rice")
dev.off()

# Let's make a bar plot of the counts just because
sp.counts <- table(flt_rpkm$Tissue.Specific)
pdf(file="Tau_Counts_AllTissues.pdf", height=8, width=8)
at <- barplot(
    sp.counts,
    ylim=c(0, max(sp.counts)*1.2),
    col="black",
    xlab="",
    ylab="Number of Genes",
    main="Counts of Tissue Specific Genes")
box()
# Put text in
text(
    x=at,
    y=max(sp.counts)*1.1,
    labels=sp.counts,
    srt=90)
dev.off()

# Let's make a heatmap of the genes above our specificity threshold.
specific <- flt_rpkm[flt_rpkm$Tau >= TAU_THRESHOLD, 1:(ncol(flt_rpkm)-3)]
specific <- as.matrix(specific)
# Put it on the log scale
specific <- log(specific+1)
# Remove genes with 0 variance
genevar <- apply(specific, 1, var)
novar <- genevar == 0
specific <- specific[!novar,]
pdf(file="WR_Tau_AllTissues_Heatmap.pdf", height=6, width=6)
pheatmap(
    specific,
    color=viridis_pal()(100),
    cluster_rows=TRUE,
    cluster_cols=TRUE,
    show_rownames=FALSE,
    treeheight_row=0,
    scale="row",
    breaks=seq(-3, 3, length.out=101))
dev.off()

# Let's try dropping root from the matrix, because we know the library is not
# very representative of root gene expression. Drop the root column, and the
# columns for Tau and gene name
flt_rpkm <- flt_rpkm[,-c(5, 10, 11)]
flt_rpkm$Tau <- as.numeric(apply(flt_rpkm, 1, calc_tau))
flt_rpkm$Gene <- as.character(raw_mat$Geneid[keep])
flt_rpkm$Tissue.Specific <- apply(flt_rpkm, 1, get_tissue)
write.csv(flt_rpkm, file="WR_Tau_WithoutRoot.csv", quote=FALSE, row.names=FALSE)

sp.counts <- table(flt_rpkm$Tissue.Specific)
pdf(file="Tau_Counts_WithoutRoot.pdf", height=8, width=8)
at <- barplot(
    sp.counts,
    ylim=c(0, max(sp.counts)*1.2),
    col="black",
    xlab="",
    ylab="Number of Genes",
    main="Counts of Tissue Specific Genes")
box()
# Put text in
text(
    x=at,
    y=max(sp.counts)*1.1,
    labels=sp.counts,
    srt=90)
dev.off()

pdf(file="Tau_Density_WithoutRoot.pdf", height=8, width=8)
tdens <- flt_rpkm$Tau[!is.na(flt_rpkm$Tau)]
hist(
    tdens,
    xlim=c(0, 1),
    breaks=100,
    col="lightgrey",
    border=NA,
    xlab="Tau",
    ylab="Count",
    main="Tissue Specificity in Northern Wild Rice\nWithout Root Sample")
dev.off()
specific <- flt_rpkm[flt_rpkm$Tau >= TAU_THRESHOLD, 1:(ncol(flt_rpkm)-3)]
specific <- as.matrix(specific)
# Put it on the log scale
specific <- log(specific+1)
# Remove genes with 0 variance
genevar <- apply(specific, 1, var)
novar <- genevar == 0
specific <- specific[!novar,]
pdf(file="WR_Tau_WithoutRoot_Heatmap.pdf", height=6, width=6)
pheatmap(
    specific,
    color=viridis_pal()(100),
    cluster_rows=TRUE,
    cluster_cols=TRUE,
    show_rownames=FALSE,
    treeheight_row=0,
    scale="row",
    breaks=seq(-3, 3, length.out=101))
dev.off()

# Maybe remove flower and sheath, too, because those are high contamination
flt_rpkm <- flt_rpkm[,-c(2, 6, 9, 10)]
flt_rpkm$Tau <- as.numeric(apply(flt_rpkm, 1, calc_tau))
flt_rpkm$Gene <- as.character(raw_mat$Geneid[keep])
flt_rpkm$Tissue.Specific <- apply(flt_rpkm, 1, get_tissue)
write.csv(flt_rpkm, file="WR_Tau_WithoutRoot_WithoutFlower_WithoutSheath.csv", quote=FALSE, row.names=FALSE)

sp.counts <- table(flt_rpkm$Tissue.Specific)
pdf(file="Tau_Counts_WithoutRoot_WithoutFlower_WithoutSheath.pdf", height=8, width=8)
at <- barplot(
    sp.counts,
    ylim=c(0, max(sp.counts)*1.2),
    col="black",
    xlab="",
    ylab="Number of Genes",
    main="Counts of Tissue Specific Genes")
box()
# Put text in
text(
    x=at,
    y=max(sp.counts)*1.1,
    labels=sp.counts,
    srt=90)
dev.off()

pdf(file="Tau_Density_WithoutRoot_WithoutFlower_WithoutSheath.pdf", height=8, width=8)
tdens <- flt_rpkm$Tau[!is.na(flt_rpkm$Tau)]
hist(
    tdens,
    xlim=c(0, 1),
    breaks=100,
    col="lightgrey",
    border=NA,
    xlab="Tau",
    ylab="Count",
    main="Tissue Specificity in Northern Wild Rice\nWithout Root, Flower, Sheath Samples")
dev.off()
specific <- flt_rpkm[flt_rpkm$Tau >= TAU_THRESHOLD, 1:(ncol(flt_rpkm)-3)]
specific <- as.matrix(specific)
# Put it on the log scale
specific <- log(specific+1)
# Remove genes with 0 variance
genevar <- apply(specific, 1, var)
novar <- genevar == 0 | is.na(genevar)
specific <- specific[!novar,]
pdf(file="WR_Tau_WithoutRoot_WithoutFlower_WithoutSheath_Heatmap.pdf", height=6, width=6)
pheatmap(
    specific,
    color=viridis_pal()(100),
    cluster_rows=TRUE,
    cluster_cols=TRUE,
    show_rownames=FALSE,
    treeheight_row=0,
    scale="row",
    breaks=seq(-3, 3, length.out=101))
dev.off()
