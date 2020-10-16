# 25 July 2020
# Starter code for reading in 4DTV results so that I can plot them

library(data.table)

#infile = sys.arg[1]

files <- Sys.glob('~/4DTV_pascal/4DTV_output/*txt')
zp_zl_files <- Sys.glob('~/4DTV_pascal/4DTV_output/zpalustris_zlatifolia/*txt')
os_os_files <- Sys.glob('~/4DTV_pascal/4DTV_output/osativa_osativa/*txt')

# Initialize an empty data.table
new_data <- data.table(NULL)

for(infile in files){
	data <- read.delim(infile, header = FALSE, sep = "\t", fill = TRUE)
	data <- as.data.table(data)
	data[, c("Col1", "Col2") := tstrsplit(V1, ":", fixed = TRUE)]
	data[,V1 := NULL]
	setnames(data, c("stat", "value"))
	data[, file := infile]
	data[, comparison := "zpalustris-zpalustris"]
	data2 <- data
	new_data <- rbind(data2[stat=="4DTV",], new_data, fill=TRUE)
}


for(infile in zp_zl_files){
	data <- read.delim(infile, header = FALSE, sep = "\t", fill = TRUE)
	data <- as.data.table(data)
	data[, c("Col1", "Col2") := tstrsplit(V1, ":", fixed = TRUE)]
	data[,V1 := NULL]
	setnames(data, c("stat", "value"))
	data[, file := infile]
	data[, comparison := "zpalustris-zlatifolia"]
	new_data <- rbind(data[stat=="4DTV",], new_data, fill=TRUE)
}

for(infile in os_os_files){
	data <- read.delim(infile, header = FALSE, sep = "\t", fill = TRUE)
	data <- as.data.table(data)
	data[, c("Col1", "Col2") := tstrsplit(V1, ":", fixed = TRUE)]
	data[,V1 := NULL]
	setnames(data, c("stat", "value"))
	data[, file := infile]
	data[, comparison := "osativa-osativa"]
	new_data <- rbind(data[stat=="4DTV",], new_data, fill=TRUE)
}


#save(new_data, file = "complete_4DTV_data.Rdata")

new_data[, value := as.numeric(value)] # Make numeric to allow rounding
new_data[, value_rounded := round(value, digits = 2)] # Round 4DTV values to 2 digits

binned <- new_data[comparison == "zpalustris-zpalustris", .N, by = "value_rounded"] # Put into bins
binned <- binned[value_rounded != 0.00] # Get rid of values of 0
binned[, as_percent := round(N/sum(binned$N), digits = 2)]
binned <- binned[order(value_rounded)] # Put in increasing order
binned <- binned[, comparison := "zpalustris-zpalustris"]

binned2 <- new_data[comparison == "zpalustris-zlatifolia", .N, by = "value_rounded"]
binned2 <- binned2[value_rounded != 0.00] # Get rid of values of 0
binned2[, as_percent := round(N/sum(binned2$N), digits = 2)]
binned2 <- binned2[order(value_rounded)] # Put in increasing order
binned2 <- binned2[, comparison := "zpalustris-zlatifolia"]

binned3 <- new_data[comparison == "osativa-osativa", .N, by = "value_rounded"]
binned3 <- binned3[value_rounded != 0.00] # Get rid of values of 0
binned3[, as_percent := round(N/sum(binned2$N), digits = 2)]
binned3 <- binned3[order(value_rounded)] # Put in increasing order
binned3 <- binned3[, comparison := "osativa-osativa"]

binned <- rbind(binned, binned2, binned3)

# Re-orer after merging
binned <- binned[order(value_rounded)] # Put in increasing order

# Set color for comparisons
binned[comparison == "zpalustris-zpalustris", col := "#007549" ]
binned[comparison == "zpalustris-zlatifolia", col := "yellow" ]
binned[comparison == "osativa-osativa", col := "black" ]

# Try out different types of plots; keep in mind some data still appear to be missing (it's not all loading)
pdf("out.pdf")
binned[comparison == "zpalustris-zpalustris", plot(value_rounded, as_percent, type="l", col=col, lwd=2, main = "4DTV", ylab= "% of paralogous pairs", xlab = "4DTV ratio", yaxt="n")]
binned[comparison == "zpalustris-zlatifolia", lines(value_rounded, as_percent, col=col, lwd=2, main = "4DTV", ylab= "% of paralogous pairs", xlab = "4DTV ratio", yaxt="n")]
binned[comparison == "osativa-osativa", lines(value_rounded, as_percent, col=col, lwd=2, main = "4DTV", ylab= "% of paralogous pairs", xlab = "4DTV ratio", yaxt="n")]
axis(2, las = 2) # makes y-axis labels easier to read
legend("topright", legend=c(expression(italic("Z. palustris - Z. palustris")), expression(italic("Z. palustris - Z. latifolia")), expression(italic("O. sativa - O. sativa"))), col=c("#007549", "yellow", "black"), lwd=3)
binned[, barplot(N, col = col, yaxt="n", ylab = "# of Occurrences", xlab = "4DTV")]
axis(2, las = 2) # makes y-axis labels easier to read)
legend("topright", legend=c(expression(italic("Z. palustris - Z. palustris")), expression(italic("Z. palustris - Z. latifolia")), expression(italic("O. sativa - O. sativa"))), col=c("#007549", "yellow", "black"), lwd=3)
dev.off()