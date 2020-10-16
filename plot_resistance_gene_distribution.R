# Load requried packages
library(graphics)
library(foreach) # required to use 'foreach()' function
library(data.table)

setwd("/Users/matthewwilliamhaas/Documents/wild_rice/disease_resistance_genes")

#  Define chromosomes and chromosome lengths (in Mb)
chromosomes = c('Chr 1', 'Chr 2', 'Chr 3', 'Chr 4', 'Chr 5', 'Chr 6', 'Chr 7', 'Chr 8',
                'Chr 9', 'Chr 10', 'Chr 11', 'Chr 12', 'Chr 13', 'Chr 14', 'Chr 15', 'Scf 16', 'Scf 458')
chr_lengths = c(95.4, 103.4, 58.8, 98.7, 66.6, 118.0, 42.6, 75.7,
                95.1, 111.4, 63.2, 105.9, 111.3, 24.0, 39.1, 13.8, 4.3)


# Define a function that adds genes to the plot
add_genes <- function(positions, xleft = 4, ybottom = NULL, xright = 6, ytop = NULL, col = NULL){
  for(j in positions){
    rect(xleft, ybottom = j+1, xright, ytop = j-1, col = col)
  }
}

# Read in data
data <- fread("updated_resistance_genes_to_plot.csv")

# Set up plots
pdf('out.pdf') # width=250, height= 250)
layout(matrix(c(1:18),2,9,byrow = TRUE))
par(mar=c(1,1,1,1))
par(oma=c(0,6,2,0))
foreach(i = chromosomes, j = chr_lengths) %do% {
    if (i == chromosomes[1] | i == chromosomes[10]){
      plot(1, type = "n", xlab = "", ylab = "", 
      xlim = c (0, 25), ylim = rev(c(0, 120)), cex = 4, axes = FALSE)
      axis(2, las = 2)
      rect(0,0,2,j, col = "grey")
      mtext(i, side = 3, at = 1)
        if (i == chromosomes[1]){
          data[chr == 1 & product_id == 1, add_genes(positions = start_mb, col='#007549')] #green
          data[chr == 1 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
          data[chr == 1 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
          data[chr == 1 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
          data[chr == 1 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
          data[chr == 1 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
          data[chr == 1 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')] 
        }
        else if (i == chromosomes[10]){
          data[chr == 10 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
          data[chr == 10 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
          data[chr == 10 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
          data[chr == 10 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
          data[chr == 10 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
          data[chr == 10 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
          data[chr == 10 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
        }
    }
  else {
      plot(1, type = "n", xlab = "", ylab = "", 
           xlim =c (0, 25), ylim = rev(c(0, 120)), axes = FALSE)
      rect(0,0,2,j, col = 'grey')
      mtext(i, side = 3, at = 1)
      
      if (i == chromosomes[2]){
        data[chr == 2 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 2 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 2 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 2 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 2 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 2 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 2 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[3]){
        data[chr == 3 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 3 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 3 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 3 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 3 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 3 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 3 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[4]){
        data[chr == 4 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 4 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 4 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 4 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 4 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 4 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 4 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[5]){
        data[chr == 5 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 5 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 5 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 5 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 5 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 5 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 5 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[6]){
        data[chr == 6 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 6 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 6 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 6 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 6 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 6 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 6 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[7]){
        data[chr == 7 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 7 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 7 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 7 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 7 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 7 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 7 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[8]){
        data[chr == 8 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 8 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 8 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 8 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 8 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 8 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 8 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[9]){
        data[chr == 9 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 9 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 9 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 9 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 9 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 9 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 9 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[10]){
        data[chr == 10 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 10 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 10 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 10 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 10 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 10 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 10 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[11]){
        data[chr == 11 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 11 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 11 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 11 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 11 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 11 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 11 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[12]){
        data[chr == 12 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 12 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 12 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 12 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 12 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 12 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 12 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
        }
      if (i == chromosomes[13]){
        data[chr == 13 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 13 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 13 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 13 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 13 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 13 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 13 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[14]){
        data[chr == 14 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 14 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 14 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 14 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 14 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 14 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 14 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[15]){
        data[chr == 15 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 15 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 15 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 15 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 15 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 15 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 15 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[16]){
        data[chr == 51 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 51 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 51 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 51 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 51 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 51 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 51 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      if (i == chromosomes[17]){
        data[chr == 453 & product_id == 1, add_genes(positions = start_mb, col='#007549')]
        data[chr == 453 & product_id == 2, add_genes(positions = start_mb, xleft = 7, xright = 9, col='red')]
        data[chr == 453 & product_id == 3, add_genes(positions = start_mb, xleft = 10, xright = 12, col='blue')]
        data[chr == 453 & product_id == 4, add_genes(positions = start_mb, xleft = 13, xright = 15, col='yellow')]
        data[chr == 453 & product_id == 5, add_genes(positions = start_mb, xleft = 16, xright = 18, col='#7130D5')] #purple
        data[chr == 453 & product_id == 6, add_genes(positions = start_mb, xleft = 19, xright = 21, col='turquoise1')]
        data[chr == 453 & product_id == 7, add_genes(positions = start_mb, xleft = 22, xright = 24, col='darksalmon')]
      }
      }
}

mtext("Physical Distance (Mb)", side = 2, outer = TRUE, line = 2) # better to add this outside of the loop
dev.off()

