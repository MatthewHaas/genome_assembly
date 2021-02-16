# Summarize the sizes of orthologous groups in Zizania palustris versus
# Oryza sativa.

ortho_sizes <- read.csv("/Users/tomkono/Dropbox/GitHub/RIS/Kimball_WR_Annotation/Results/Orthofinder/Zp_vs_Os_Orthologue_Count.csv", header=TRUE)

# Define a function that returns a label that classifies the orthologous group
# by gene count in both species
classify_og <- function(x) {
    zp <- as.numeric(x[2])
    os <- as.numeric(x[3])
    if(zp == os) {
        return("Equal")
    } else if(zp > os) {
        return("Expanded in Z. palustris")
    } else if(zp < os) {
        return("Expanded in O. sativa")
    } else {
        return(NA)
    }
}

ortho_class <- apply(ortho_sizes, 1, classify_og)

# Plot the sizes!
pdf(file="Zp_vs_Os_Orthologue_Counts.pdf", height=6, width=6)
plot(
    jitter(ortho_sizes$Oryza_sativa) ~ (ortho_sizes$Zizania_palustris),
    cex=0.5,
    xlab="Orthologus Genes in Northern Wild Rice",
    ylab="Orthologus Genes in Asian Rice",
    main="")
abline(a=0, b=1)
dev.off()

pdf(file="Zp_vs_Os_Expansions.pdf", height=6, width=4)
par(mar=c(7.5, 3, 1, 1))
cts <- table(ortho_class)
names(cts) <- NULL
at <- barplot(cts, col=c("black", "#1f78b4", "#b2df8a"), axes=FALSE, ylim=c(0, 10000))
axis(side=2, cex.axis=0.75)
axis(
    side=1,
    at=at,
    labels=c("Equal", "O. sativa Expanded", "Z. palustris Expanded"),
    las=2,
    cex.axis=0.75)
box()
dev.off()
