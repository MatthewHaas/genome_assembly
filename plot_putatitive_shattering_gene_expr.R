# 27 July 2020
# Updated: 5 October 2020
# Make plots for shattering genes

# Load required packages
library(data.table)
library(reshape2)

infile = args[1] # current file = cpm_list.txt

# Read in data
x <- fread(infile)

y <- melt(x, id.vars = 1, variable.name = "tissue")

# We think _sh4_ is FUN028835 AND _sh8_ is FUN_004935
# With the new gene names, this corresponds to ZPchr0458G22499 (FUN_028835) and ZPchr0008g13714 (FUN_004935)
# But we are not sure about the other shattering genes (_qsh1_, _sh3_, and _sh5_).
#putative_shattering_genes <- c("FUN_022397", "FUN_010896", "FUN_028711", "FUN_021597", "FUN_021815",
#							   "FUN_019929", "FUN_004935", "FUN_028835")

putative_shattering_genes <- c("FUN_028835", "FUN_023859", "FUN_038639", "FUN_038661", "FUN_011199", 
                               "FUN_015882", "FUN_014124", "FUN_038639")
# "FUN_013657" was removed; not in data set
new_names <- c("ZPchr0458g22499", "ZPchr0003g18426", "ZPchr0010g10516", "ZPchr0010g7757", "ZPchr0013g34051",
               "ZPchr0001g31104", "ZPchr0005g15825", "ZPchr0010g10516")
# "ZPchr0005g1585" removed because it corresponds to FUN_013657


# FUN_028835 (sh4), FUN_023859 (Sh1), FUN_038661 (qSH1), FUN_011199 (SHAT1), FUN_038661 (sh5)
# others: FUN_038639, FUN_015882, FUN_013657, FUN_038639
# ZPchr0458g22499 (sh4), ZPchr0003g18426 (Sh1), ZPchr0010g7757 (qSH1), ZPchr0013g34051 (SHAT1), ZPchr0010g7757 (sh5)
# others: ZPchr0010g10516, ZPchr0001g31104, ZPchr0005g1585, ZPchr0010g10516

#new_names <- c("ZPchr0003g18492", "ZPchr0013g36816", "ZPchr0458g22566", "ZPchr0004g38673",
#              "ZPchr0004g39424", "ZPchr0004g38210", "ZPchr0008g13714", "ZPchr0458g22499")

# Make the plot
# The variable j allows for the translation of the old gene names (FUN*) with the new names (ZP*).
# I still look up by the old FUN names because that's how they are given in the cpm_list.txt file,
# but we have a translation table for both genes and chromosomes, so the information is interchangeable.
pdf("201014_putative_shattering_genes.pdf")
layout(matrix(c(1:8),2,4,byrow = TRUE))
j = 1
for(i in putative_shattering_genes){
	counts <- y[genes == i,]$value
	par(cex.main = 0.9)
	par(family = 'serif') # Makes font Times New Roman-like (serif)
	if (i == "FUN_028835"){
	  barplot(counts, main = expression(atop("ZPchr0458g22499", italic('sh4'))), font.main = 2, ylab = "Relative expression", ylim = c(-6.5, 6.5), 
	          names.arg = unique(y$tissue), las = 2, col = "steelblue")
  j = j + 1
	}
	else if (i == "FUN_023859"){
	  barplot(counts, main = expression(atop("ZPchr0003g18426", italic('Sh1'))), ylab = "Relative expression", ylim = c(-6.5, 6.5), 
	          names.arg = unique(y$tissue), las = 2, col = "steelblue")
	  j = j + 1
	}
	else if (i == "FUN_038661"){
	    barplot(counts, main = expression(atop("ZPchr0010g7757", italic('qSH1/sh5'))), ylab = "Relative expression", ylim = c(-6.5, 6.5), 
	            names.arg = unique(y$tissue), las = 2, col = "steelblue")
	    j = j + 1
	}
	else if (i == "FUN_011199"){
	    barplot(counts, main = expression(atop("ZPchr0013g34051",paste("(", italic('SHAT1'), ")"))), ylab = "Relative expression", ylim = c(-6.5, 6.5), 
	            names.arg = unique(y$tissue), las = 2, col = "steelblue")
	    j = j + 1
	}
	else {
	  print(i)
	  barplot(counts, main = new_names[j], ylab = "Relative expression", ylim = c(-6.5, 6.5), 
	          names.arg = unique(y$tissue), font.main = 1, las = 2, col = "steelblue")
	  j = j + 1
	}
}
dev.off()


pdf("201014_putative_shattering_genes.pdf")
layout(matrix(c(1:8),2,4,byrow = TRUE))
j = 1
for(i in putative_shattering_genes){
  counts <- y[genes == i,]$value
  par(cex.main = 0.9)
  par(family = 'serif') # Makes font Times New Roman-like (serif)
  if (i == "FUN_028835"){
    barplot(counts, main = expression(atop("ZPchr0458g22499", italic('sh4'))), font.main = 2, ylab = "Relative expression", ylim = c(-6.5, 6.5), 
            names.arg = unique(y$tissue), las = 2, col = "steelblue")
    j = j + 1
  }
  else if (i == "FUN_023859"){
    barplot(counts, main = expression(atop("ZPchr0003g18426", italic('Sh1'))), ylab = "Relative expression", ylim = c(-6.5, 6.5), 
            names.arg = unique(y$tissue), las = 2, col = "steelblue")
    j = j + 1
  }
  else if (i == "FUN_038661"){
    barplot(counts, main = expression(atop("ZPchr0010g7757", italic('qSH1/sh5'))), ylab = "Relative expression", ylim = c(-6.5, 6.5), 
            names.arg = unique(y$tissue), las = 2, col = "steelblue")
    j = j + 1
  }
  else if (i == "FUN_011199"){
    barplot(counts, main = expression(atop("ZPchr0013g34051",paste("(", italic('SHAT1'), ")"))), ylab = "Relative expression", ylim = c(-6.5, 6.5), 
            names.arg = unique(y$tissue), las = 2, col = "steelblue")
    j = j + 1
  }
  else {
    print(i)
    barplot(counts, main = new_names[j], ylab = "Relative expression", ylim = c(-6.5, 6.5), 
            names.arg = unique(y$tissue), font.main = 1, las = 2, col = "steelblue")
    j = j + 1
  }
}
dev.off()

