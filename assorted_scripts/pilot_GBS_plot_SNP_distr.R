library(data.table)

setwd("~/Documents/wild_rice")

x <- fread("210302_samtools_ZPchrs_filtered_concat.csv")

scaffolds <- unique(x$scaffold)


pdf("210302_snp_distribution.pdf")
layout(matrix(c(1,5,9,13,17,
                2,6,10,14,0,
                3,7,11,15,0,
                4,8,12,16,0), nrow=5), widths=c(1,1,1,1))
par(oma=c(3,6,3,0))
par(mar=c(2,2,2,2))
for (i in scaffolds){
  y <- x[scaffold == i, .N, key=.(bin=position %/% 1e6)]
  y[, plot(x = bin, y = N, main = i,
           xlab = "chromosome position (Mb)",
           ylab = "# of SNPs",
           family = "serif",
           type = "l",
           las = 1)]
 mtext(text = "Number of SNPs", family = "serif", side = 2, line = 1, outer = TRUE)  
 mtext(text = "Chromosome position (Mb)", family = "serif", side = 1, outer = TRUE)
}
dev.off()
