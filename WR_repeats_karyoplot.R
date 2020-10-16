##----- karyoplotR - create BSgenome object-----------##
# make zizania palustris seed file - manually
# add seqnames to the seed file with the make.vector.py script. Produces chrom_vector.txt file to append to the end of the seed file.
# cat zizania_palustris-seed chrom_vector.txt > zizania_palustris-seed2
setwd("/home/jkimball/shared/WR_Annotation/genome/zizania_palustris_scaffolds")
library(BSgenome)
Zpalus_seed <- "/home/jkimball/shared/WR_Annotation/genome/zizania_palustris-seed2"
forgeBSgenomeDataPkg(Zpalus_seed)

# exit R
cd /home/jkimball/shared/WR_Annotation/genome/zizania_palustris_scaffolds/
R CMD build BSgenome.Zpalus/
# for some reason, the check function threw an error saying that the disk quota has been exceeded. No idea why that is.
R CMD check BSgenome.Zpalus_1.0.0.tar.gz
R CMD INSTALL BSgenome.Zpalus_1.0.0.tar.gz 

##----- karyoplotR ----------------------------------##
# module load R/3.6.0
library(BSgenome.Zpalus)
library(karyoploteR)
library(scales)
library(GenomicFeatures)
Zpalus <- BSgenome.Zpalus

# have to rename the seq levels to be either ENSEMBL, UCSC -- some recognizable format
#Zpalus2 <- renameSeqlevels(Zpalus,sub("Scaffold_","chr",seqlevels(Zpalus)))
# This is stupid and messy, but without renaming the Zpalus object each time, the object is overwritten
Zpalus2 <- renameSeqlevels(Zpalus,sub("^Scaffold_1$", "chr9",seqlevels(Zpalus)))
Zpalus3 <- renameSeqlevels(Zpalus,sub("^Scaffold_3$", "chr3",seqlevels(Zpalus2)))
Zpalus4 <- renameSeqlevels(Zpalus,sub("^Scaffold_7$", "chr15",seqlevels(Zpalus3)))
Zpalus5 <- renameSeqlevels(Zpalus,sub("^Scaffold_9$", "chr11",seqlevels(Zpalus4)))
Zpalus6 <- renameSeqlevels(Zpalus,sub("^Scaffold_13$", "chr1",seqlevels(Zpalus5)))
Zpalus7 <- renameSeqlevels(Zpalus,sub("^Scaffold_18$", "chr4",seqlevels(Zpalus6)))
Zpalus8 <- renameSeqlevels(Zpalus,sub("^Scaffold_48$", "chr6",seqlevels(Zpalus7)))
Zpalus9 <- renameSeqlevels(Zpalus,sub("^Scaffold_51$", "chr51",seqlevels(Zpalus8)))
Zpalus10 <- renameSeqlevels(Zpalus,sub("^Scaffold_70$", "chr10",seqlevels(Zpalus9)))
Zpalus11 <- renameSeqlevels(Zpalus,sub("^Scaffold_93$", "chr2",seqlevels(Zpalus10)))
Zpalus12 <- renameSeqlevels(Zpalus,sub("^Scaffold_415$", "chr12",seqlevels(Zpalus11)))
Zpalus13 <- renameSeqlevels(Zpalus,sub("^Scaffold_453$", "chr453",seqlevels(Zpalus12)))
Zpalus14 <- renameSeqlevels(Zpalus,sub("^Scaffold_693$", "chr14",seqlevels(Zpalus13)))
Zpalus15 <- renameSeqlevels(Zpalus,sub("^Scaffold_1062$", "chr8",seqlevels(Zpalus14)))
Zpalus16 <- renameSeqlevels(Zpalus,sub("^Scaffold_1063$", "chr7",seqlevels(Zpalus15)))
Zpalus17 <- renameSeqlevels(Zpalus,sub("^Scaffold_1064$", "chr13",seqlevels(Zpalus16)))
Zpalus18 <- renameSeqlevels(Zpalus,sub("^Scaffold_1065$", "chr5",seqlevels(Zpalus17)))

# Make object Zpalus2 again
Zpalus2 <- Zpalus18

seqlevels(Zpalus2)
seqlevelsStyle(Zpalus2) <- "UCSC"
kp <- plotKaryotype(genome=Zpalus2)

# number of WR scaffolds > 500,000
sum(seqlengths(Zpalus2) > 500000)
# 17

top17_scaffolds <- seqnames(Zpalus2)[seqlengths(Zpalus2) > 500000]

top17_chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
					   "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr51", "chr453")

#kp <- plotKaryotype(genome=Zpalus2,chromosomes=top17_scaffolds)
# I think this needs to be altered in order to work.. rename top 17 scaffolds
kp <- plotKaryotype(genome=Zpalus2,chromosomes=top17_chromosomes)

# load genes
txdb <- makeTxDbFromGFF("/home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/predict_results/rice.gff3", 
                        format = "gff3",
                        dataSource="none",
                        organism="Zizania palustris")
all.genes <- genes(txdb)
# change the seqlevels to match the BSgenome object
#seqlevels(all.genes) <- sub("Scaffold_","chr",seqnames(all.genes)) #should seqnames() be seqlevels()???

# Trying to change the individual scaffold names with actual chromosome names like I did with repeats
# First two seem to work; returns error for chr 15.
seqlevels(all.genes)  <- sub("^Scaffold_1$", "chr9", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_3$", "chr3", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_7$", "chr15", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_9$", "chr11", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_13$", "chr1", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_18$", "chr4", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_48$", "chr6", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_51$", "chr51", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_70$", "chr10", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_93$", "chr2", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_415$", "chr12", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_453$", "chr453", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_693$", "chr14", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_1062$", "chr8", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_1063$", "chr7", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_1064$", "chr13", seqlevels(all.genes))
seqlevels(all.genes) <- sub("^Scaffold_1065$", "chr5", seqlevels(all.genes))

# load repeats
repeats <- read.delim("/home/jkimball/shared/WR_Annotation/genome/zizania_palustris_13Nov2018_okGsv.bed",header=FALSE,sep="\t")
colnames(repeats) <- c("seqid","start","end","repeats")

# rename scaffolds to match seqids from the BS genome object
#repeats$seqid <- sub("Scaffold_","chr",repeats$seqid)
	
repeats$seqid <- sub("^chr1$", "chr9", repeats$seqid)
repeats$seqid <- sub("^chr3$", "chr3", repeats$seqid)
repeats$seqid <- sub("^chr7$", "chr15", repeats$seqid)
repeats$seqid <- sub("^chr9$", "chr11", repeats$seqid)
repeats$seqid <- sub("^chr13$", "chr1", repeats$seqid)
repeats$seqid <- sub("^chr18$", "chr4", repeats$seqid)
repeats$seqid <- sub("^chr48$", "chr6", repeats$seqid)
repeats$seqid <- sub("^chr51$", "chr51", repeats$seqid)
repeats$seqid <- sub("^chr70$", "chr10", repeats$seqid)
repeats$seqid <- sub("^chr93$", "chr2", repeats$seqid)
repeats$seqid <- sub("^chr415$", "chr12", repeats$seqid)
repeats$seqid <- sub("^chr453$", "chr453", repeats$seqid)
repeats$seqid <- sub("^chr693$", "chr14", repeats$seqid)
repeats$seqid <- sub("^chr1062$", "chr8", repeats$seqid)
repeats$seqid <- sub("^chr1063$", "chr7", repeats$seqid)
repeats$seqid <- sub("^chr1064$", "chr13", repeats$seqid)
repeats$seqid <- sub("^chr1065$", "chr5", repeats$seqid)
	
# subset repeats - LTR, LINE, SINE, and DNA
LTR <- repeats[repeats$repeats == "LTR",1:3]
LINE <- repeats[repeats$repeats == "LINE",1:3]
SINE <- repeats[repeats$repeats == "SINE",1:3]
DNA <- repeats[repeats$repeats == "DNA",1:3]

# make repeat file into GR ranges object
LTR <- makeGRangesFromDataFrame(LTR)
LINE <- makeGRangesFromDataFrame(LINE)
SINE <- makeGRangesFromDataFrame(SINE)
DNA <- makeGRangesFromDataFrame(DNA)
allreps <- makeGRangesFromDataFrame(repeats)

# Change chromosome names for all repetitive elements
# LTRs
seqlevels(LTR)  <- sub("^Scaffold_1$", "chr9", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_3$", "chr3", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_7$", "chr15", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_9$", "chr11", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_13$", "chr1", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_18$", "chr4", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_48$", "chr6", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_51$", "chr51", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_70$", "chr10", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_93$", "chr2", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_415$", "chr12", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_453$", "chr453", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_693$", "chr14", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_1062$", "chr8", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_1063$", "chr7", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_1064$", "chr13", seqlevels(LTR))
seqlevels(LTR) <- sub("^Scaffold_1065$", "chr5", seqlevels(LTR))
# LINEs
seqlevels(LINE)  <- sub("^Scaffold_1$", "chr9", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_3$", "chr3", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_7$", "chr15", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_9$", "chr11", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_13$", "chr1", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_18$", "chr4", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_48$", "chr6", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_51$", "chr51", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_70$", "chr10", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_93$", "chr2", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_415$", "chr12", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_453$", "chr453", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_693$", "chr14", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_1062$", "chr8", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_1063$", "chr7", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_1064$", "chr13", seqlevels(LINE))
seqlevels(LINE) <- sub("^Scaffold_1065$", "chr5", seqlevels(LINE))
# SINEs
seqlevels(SINE)  <- sub("^Scaffold_1$", "chr9", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_3$", "chr3", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_7$", "chr15", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_9$", "chr11", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_13$", "chr1", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_18$", "chr4", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_48$", "chr6", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_51$", "chr51", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_70$", "chr10", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_93$", "chr2", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_415$", "chr12", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_453$", "chr453", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_693$", "chr14", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_1062$", "chr8", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_1063$", "chr7", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_1064$", "chr13", seqlevels(SINE))
seqlevels(SINE) <- sub("^Scaffold_1065$", "chr5", seqlevels(SINE))
# DNA
seqlevels(DNA)  <- sub("^Scaffold_1$", "chr9", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_3$", "chr3", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_7$", "chr15", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_9$", "chr11", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_13$", "chr1", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_18$", "chr4", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_48$", "chr6", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_51$", "chr51", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_70$", "chr10", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_93$", "chr2", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_415$", "chr12", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_453$", "chr453", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_693$", "chr14", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_1062$", "chr8", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_1063$", "chr7", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_1064$", "chr13", seqlevels(DNA))
seqlevels(DNA) <- sub("^Scaffold_1065$", "chr5", seqlevels(DNA))

# How to reorder?					   
seqlevels(Zpalus2) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
					   "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr51", "chr453")

kp <- plotKaryotype(genome=Zpalus2,chromosomes=top17_chromosomes)
kpPlotDensity(kp, LTR,col="blue")
kpAxis(kp)
#kpPoints(kp, data=LINE) #throws error. no Y value??
kpPlotDensity(kp, LINE,col = alpha("grey", 0.4))
kpPlotDensity(kp, LINE,col = "grey")
kpPlotDensity(kp, DNA,col = "red")
kpPlotDensity(kp, SINE,col = "yellow")
kpPlotDensity(kp,all.genes,col="purple")
#kpPlotDensity(kp,allreps,col="green")

