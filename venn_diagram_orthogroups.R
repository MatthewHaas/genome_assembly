# 21 December 2020
# Purpose of this code is to make a venn diagram for overlapping Orthogroups for the NWR genome paper
# From OrthoFinder run with only the species featured in the venn diagram

setwd("~/Documents/wild_rice")

# Load required packages
library(data.table)
library(VennDiagram)
# Since we wanted to add commas as the separator for thousands, I modified the venn.diagram function.
# We need to use source() to load the modified function.
source('~/Documents/wild_rice/scripts/modified_venn_diagram.R')

# Read in data
x <- fread("tab-separated_value_files/Orthogroups_GeneCount_fig2.tsv") #downloaded to my computer and name changed

# Find all rows for which each species has at least one gene present in the orthogroup
O_sativa <- x[Oryza_sativa > 0]
Z_palustris <- x[Zizania_palustris > 0]
Z_mays <- x[Zea_mays > 0]
B_distachyon <- x[Brachypodium_distachyon > 0]
S_bicolor <- x[Sorghum_bicolor > 0]

# Get the orthogroup IDS for each species.
Osativa_orthogroups <- O_sativa$Orthogroup
Zpalustris_orthogroups <- Z_palustris$Orthogroup
Zmays_orthogroups <- Z_mays$Orthogroup
Bdistachyon_orthogroups <- B_distachyon$Orthogroup
Sbicolor_orthogroups <- S_bicolor$Orthogroup

# Save data
save(Osativa_orthogroups, Zpalustris_orthogroups, Lperrieri_orthogroups, Zmays_orthogroups, Bdistachyon_orthogroups, Sbicolor_orthogroups, file="orthogroup_names_for_venn_diagram.Rdata")

# In RStudio now... I could not get the venn.diagram function to work. Could not write to PNG for some reason.
#setwd("~/Documents/wild_rice/tab-separated_value_files")
#load("orthogroup_names_for_venn_diagram.Rdata") # File contains sets of Orthogroups to feed into the venn.diagram function.

colors= c("#6b7fff", "#c3db0f", "#ff4059", "#2cff21", "#de4dff")

# Make the Venn Diagram
venn.diagram(x=list(Osativa_orthogroups, Zpalustris_orthogroups, Zmays_orthogroups, Bdistachyon_orthogroups, Sbicolor_orthogroups),
			 category.names = c(expression(atop(italic("O. sativa"), plain("19,526"))), expression(atop(italic("Z. palustris"), plain("18,589"))),
			 expression(atop(italic("Z. mays"), plain("24,978"))), expression(atop(italic("B. distachyon"),plain("19,142"))),
			 expression(atop(italic("S. bicolor"), plain("19,849")))),
			 filename = "venn_diagram_figure_2.png",
			 output=TRUE,
			 imagetype="png",
			 scaled = FALSE,
			 col = "black",
			 fill = colors,
			 cat.col = colors,
			 #cat.just = list(c(0.5, 0), c(1,1), c(1,1), c(1,1), c(1,1)),
			 cat.cex = 0.95,
			 cat.dist = 0.35,
	     		 margin = 0.30)
