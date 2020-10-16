# 23 June 2020
# Purpose of this code is to make a venn diagram for overlapping Orthogroups for the NWR genome paper
# This script uses the updated OrthoFinder results that includes Z. latifolia and other Oryza species.
# Path to data: /home/jkimball/shared/WR_Annotation/Orthofinder_2020-06/Grass_Proteomes/OrthoFinder/Results_Jun02/WorkingDirectory/OrthoFinder/Results_Jun18
# Done in RStudio

# Set working directory
setwd("~/Documents/wild_rice/tab-separated_value_files")

# Load required packages
library(data.table)
library(VennDiagram)
# Since we wanted to add commas as the separator for thousands, I modified the venn.diagram function.
# We need to use source() to load the modified function.
source('~/Documents/wild_rice/scripts/modified_venn_diagram.R')

# Read in data
x <- fread("Orthogroups.GeneCount.tsv")

# Find all rows for which each species has at least one gene present in the orthogroup
Z_palustris <- x[Zizania_palustris > 0]
Z_latifolia <- x[Zizania_latifolia > 0]
O_sativa <- x[Oryza_sativa > 0]
O_glaberrima<- x[Oryza_glaberrima > 0]
O_nivara <- x[Oryza_nivara > 0]
O_rufipogon <- x[Oryza_rufipogon > 0]
O_barthii <- x[Oryza_barthii > 0]

# Get the orthogroup IDS for each species.
Zpalustris_orthogroups <- Z_palustris$Orthogroup
Zlatifolia_orthogroups <- Z_latifolia$Orthogroup
Osativa_orthogroups <- O_sativa$Orthogroup
Oglaberrima_orthogroups <- O_glaberrima$Orthogroup
Onivara_orthogroups <- O_nivara$Orthogroup
Orufipogon_orthogroups <- O_rufipogon$Orthogroup
Obarthii_orthogroups <- O_barthii$Orthogroup

# Save data
save(Zpalustris_orthogroups, Zlatifolia_orthogroups, Osativa_orthogroups,
Oglaberrima_orthogroups, Onivara_orthogroups, Orufipogon_orthogroups,
Obarthii_orthogroups, file="orthogroup_names_for_venn_diagram_update_June_2020.Rdata")

colors= c("#6b7fff", "#c3db0f", "#ff4059", "#2cff21", "#de4dff")

# Make the Venn Diagram
venn.diagram(x=list(Osativa_orthogroups, Zpalustris_orthogroups, Zlatifolia_orthogroups, Oglaberrima_orthogroups, Orufipogon_orthogroups),
			 category.names = c(expression(atop(italic("O. sativa"), plain("20,309"))), expression(atop(italic("Z. palustris"), plain("17,536"))),
			 expression(atop(italic("Z. latifolia"), plain("17,168"))), expression(atop(italic("O. glaberrima"),plain("19,653"))),
			 expression(atop(italic("O. rufipogon"), plain("23,880")))),
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
