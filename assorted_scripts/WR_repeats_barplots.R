###------Kimball WR repeat content ---------###

library(RColorBrewer)
pal = brewer.pal(n = 12, name = "Paired")
setwd("/home/jkimball/shared/WR_Annotation/genome")
d <- read.delim("zizania_palustris_13Nov2018_okGsv.fasta.tab.out",row.names=NULL,header=FALSE)
# repeat length
d$rlen <- d$V7-d$V6

# Top 20 abundant repeats
par(mar=c(10,6,2,2))
barplot(sort(table(d$V11),decreasing=T)[1:20],las=2)


# Cumulative repeats lengths
test <- d[,c(11,17)]
test2 <- aggregate(test$rlen,by=list(test$V11),sum)

# pie chart of cumulative length
par(mar=c(4,4,4,8))
pie(test2$x,labels=test2$Group.1,col=pal)

# bar chart of the cumulative repeat length
test2 <- test2[order(test2$x,decreasing=T),]
par(mar=c(8,8,4,4))
barplot(test2$x[1:20],las=2,names.arg=test2$Group.1[1:20])

#barchart showing average length of repeats
test <- d[,c(11,17)]
#test3 <- aggregate(test$rlen,by=list(test$V11),mean)
#test3 <- test3[order(test3$x,decreasing=T),]
#barplot(test3$x,las=2,names.arg=test2$Group.1)
boxplot(rlen~V11,data=test,las=2)

# barplot showing % of genome
barplot((test2$x[1:25]/1288768912)*100,las=2,names.arg=test2$Group.1[1:25],ylab="% of genome",ylim=c(0,40))
