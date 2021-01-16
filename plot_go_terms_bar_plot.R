# 10 July 2020

# set working directory
setwd("/Users/matthewwilliamhaas/Documents/wild_rice/annotation_figures")

# load required pacakges
library(data.table)

# read in data
cc <- fread("go_cellular_component_genes.csv")
mf <- fread("go_molecular_function_genes.csv")
bp <- fread("go_biological_process_genes.csv")

# change column names
cc <- setnames(cc, c("gene", "annotation", 'GO_term'))
mf <- setnames(mf, c("gene", "annotation", 'GO_term'))
bp <- setnames(bp, c("gene", "annotation", 'GO_term'))

# replace extraneous characters from GO_terms column
cc <- cc[, GO_term := sub("C:", "", GO_term)]
mf <- mf[, GO_term := sub("F:", "", GO_term)]
bp <- bp[, GO_term := sub("P:", "", GO_term)]

cc <- cc[, GO_term := sub(";", "", GO_term)]
mf <- mf[, GO_term := sub(";", "", GO_term)]
bp <- bp[, GO_term := sub(";", "", GO_term)]

# count occurrences of each GO term
cc <- cc[, .(.N), by="GO_term"]
mf <- mf[, .(.N), by="GO_term"]
bp <- bp[, .(.N), by="GO_term"]

# put in decreasing order
cc <- cc[order(-N)]
mf <- mf[order(-N)]
bp <- bp[order(-N)]

# define top 13 GO terms
top13_cc_go_terms <- cc[c(1:13)]$GO_term
top13_mf_go_terms <- mf[c(1:13)]$GO_term
top13_bp_go_terms <- bp[c(1:13)]$GO_term


# make bar plot and save to file
png("out.png")
layout(matrix(c(1,2,3), nrow=1))
cc[GO_term %in% top13_cc_go_terms, barplot(cc$N[c(1:13)], main="Cellular Component", yaxt='n')]
axis(2, las=2)
mf[GO_term %in% top13_mf_go_terms, barplot(mf$N[c(1:13)], main="Molecular Function", yaxt='n')]
axis(2, las=2)
bp[GO_term %in% top13_bp_go_terms, barplot(bp$N[c(1:13)], main="Biological Process", yaxt='n')]
axis(2, las=2)
dev.off()