##----- karyoplotR - create BSgenome object-----------##
# make zizania palustris seed file - manually
# add seqnames to the seed file with the make.vector.py script. Produces chrom_vector.txt file to append to the end of the seed file.
# cat zizania_palustris-seed chrom_vector.txt > zizania_palustris-seed2
setwd("/home/jkimball/shared/WR_Annotation/genome/zizania_palustris_scaffolds")
library(BSgenome)
Zpalus_seed <- "/home/jkimball/shared/WR_Annotation/genome/zizania_palustris-seed2"
forgeBSgenomeDataPkg(Zpalus_seed)

# exit R
# cd /home/jkimball/shared/WR_Annotation/genome/zizania_palustris_scaffolds/
# R CMD build BSgenome.Zpalus/
# for some reason, the check function threw an error saying that the disk quota has been exceeded. No idea why that is.
# R CMD check BSgenome.Zpalus_1.0.0.tar.gz
# R CMD INSTALL BSgenome.Zpalus_1.0.0.tar.gz 

##----- karyoplotR ----------------------------------##
# module load R/3.6.0
library(BSgenome.Zpalus)
library(karyoploteR)
library(scales)
Zpalus <- BSgenome.Zpalus

# Rename the seq levels to be either ENSEMBL, UCSC -- some recognizable format
palus2 <- renameSeqlevels(Zpalus,sub("Scaffold_","chr",seqlevels(Zpalus)))
seqlevels(Zpalus2)
seqlevelsStyle(Zpalus2) <- "UCSC"
kp <- plotKaryotype(genome=Zpalus2)

# number of WR scaffolds > 500,000
sum(seqlengths(Zpalus2) > 500000)
# 17

top17_scaffolds <- seqnames(Zpalus2)[seqlengths(Zpalus2) > 500000]

kp <- plotKaryotype(genome=Zpalus2,chromosomes=top17_scaffolds)
# load genes
txdb <- makeTxDbFromGFF("/home/jkimball/shared/WR_Annotation/Annotation/wild_rice3/predict_results/rice.gff3", 
                        format = "gff3",
                        dataSource="none",
                        organism="Zizania palustris")
all.genes <- genes(txdb)
# change the seqlevels to match the BSgenome object
seqlevels(all.genes) <- sub("Scaffold_","chr",seqnames(all.genes))

# load repeats
repeats <- read.delim("/home/jkimball/shared/WR_Annotation/genome/zizania_palustris_13Nov2018_okGsv.bed",header=FALSE,sep="\t")
colnames(repeats) <- c("seqid","start","end","repeats")

# rename scaffolds to match seqids from the BS genome object
repeats$seqid <- sub("Scaffold_","chr",repeats$seqid)

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

kp <- plotKaryotype(genome=Zpalus2,chromosomes=top17_scaffolds)
kpPlotDensity(kp,LTR,col="blue")
kpAxis(kp)
kpPoints(kp, data=LINE)
kpPlotDensity(kp, LINE,col = alpha("gray", 0.4))
kpPlotDensity(kp, LINE,col = "gray")
kpPlotDensity(kp, DNA,col = "red")
kpPlotDensity(kp, SINE,col = "yellow")
kpPlotDensity(kp,all.genes,col="purple")
kpPlotDensity(kp,allreps,col="green")
