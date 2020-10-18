# README for genome_assembly
## Annotation of the Zizania palustris genome


## dotplot.py
This script was originally written by Haibao Tang (https://github.com/tanghaibao/jcvi). I am including the script here because I modified it in order to create my plots. The following changes were made by hard-coding my desired output into the original script: 1) The font color of the chromosome labels and positions were changed from grey to black, 2) the labels for the x and y axes were changed to _Zizania palustris_ and _Oryza sativa_ (respectively) rather than 'wild_rice' and 'oryza' (which are the BED file names), and 3) the xlimit was slightly increased (along with the length of chr 15, scf 16, and scf 458 (in order to make the chromosome labels legible).

## karyotpe.py
This script was originally written by Haibao Tang (https://github.com/tanghaibao/jcvi). I am including the script here because I modified it in order to create my plots. Line 40 was altered so that **arg[5]** (the name we assign to each track in the layout file) is printed in italics. Line 239 was also changed (dividing vpad by 2 was removed) to make extra room on the margin so that _Zizania palustris_ could be fully written out (versus abbreviating it as _Z. palustris_--which also didn't fit initally--it ran into the representations of the chromosomes.
