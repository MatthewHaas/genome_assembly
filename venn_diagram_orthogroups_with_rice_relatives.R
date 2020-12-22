# 21 December 2020
# Purpose of this code is to make a venn diagram for overlapping Orthogroups for the NWR genome paper
# This script uses the updated OrthoFinder results that includes Z. latifolia and other Oryza species.
# Done in RStudio

# Set working directory
setwd("~/Documents/wild_rice")

# Load required packages
library(data.table)
library(VennDiagram)
# Since we wanted to add commas as the separator for thousands, I modified the venn.diagram function.
# We need to use source() to load the modified function.
source('~/Documents/wild_rice/scripts/modified_venn_diagram.R')

# Read in data
x <- fread("tab-separated_value_files/Orthogroups_GeneCount_rice_relatives.tsv")

# Find all rows for which each species has at least one gene present in the orthogroup
Z_palustris <- x[Zizania_palustris > 0]
Z_latifolia <- x[Zizania_latifolia > 0]
O_sativa <- x[Oryza_sativa > 0]
O_glaberrima<- x[Oryza_glaberrima > 0]
O_rufipogon <- x[Oryza_rufipogon > 0]

# Get the orthogroup IDS for each species.
Zpalustris_orthogroups <- Z_palustris$Orthogroup
Zlatifolia_orthogroups <- Z_latifolia$Orthogroup
Osativa_orthogroups <- O_sativa$Orthogroup
Oglaberrima_orthogroups <- O_glaberrima$Orthogroup
#Onivara_orthogroups <- O_nivara$Orthogroup
Orufipogon_orthogroups <- O_rufipogon$Orthogroup
#Obarthii_orthogroups <- O_barthii$Orthogroup

# Save data
#save(Zpalustris_orthogroups, Zlatifolia_orthogroups, Osativa_orthogroups,
#Oglaberrima_orthogroups, Onivara_orthogroups, Orufipogon_orthogroups,
#Obarthii_orthogroups, file="orthogroup_names_for_venn_diagram_update_June_2020.Rdata")

colors= c("#6b7fff", "#c3db0f", "#ff4059", "#2cff21", "#de4dff")

# Make the Venn Diagram
venn.diagram(x=list(Osativa_orthogroups, Zpalustris_orthogroups, Zlatifolia_orthogroups, Oglaberrima_orthogroups, Orufipogon_orthogroups),
			 category.names = c(expression(atop(italic("O. sativa"), plain("24,940"))), expression(atop(italic("Z. palustris"), plain("20,899"))),
			 expression(atop(italic("Z. latifolia"), plain("20,457"))), expression(atop(italic("O. glaberrima"),plain("24,328"))),
			 expression(atop(italic("O. rufipogon"), plain("27,136")))),
			 filename = "venn_diagram_with_rice_relatives.png",
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
