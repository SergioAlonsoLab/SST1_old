library(gdata)
library(magrittr)
library(ggplot2)


# bluepalette

bp <- colorRampPalette(c("darkblue","lightblue"))

cellcycle <- read.xls("data/CellCycleBEA.xlsx")

d0 <- t(apply(cellcycle[,2:4],1,function(x) x/sum(x))) * 100
rownames(d0) <- cellcycle$Sample
d0 <- as.data.frame(d0)
d0 <- d0[c(15,1:14),]
colnames(d0) <- c("G1","S","G2")
par(mar=c(5,8,2,2))
barplot(t(d0[15:1,1:3]),
        ylim=c(0.5,20),
        horiz = T,las=1,
        col=bp(3)) -> mids

text(d0$G1/2,rev(mids),sprintf("%1.0f%%",d0$G1),col="white")
text(d0$G1+d0$S/2,rev(mids),sprintf("%1.0f%%",d0$S),col="white")
text(d0$G1+d0$S+d0$G2/2,rev(mids),sprintf("%1.0f%%",d0$G2),col="black")
title(xlab="Percentage of cells")
legend("topright",legend = c("G1","S","G2"),fill = bp(3),horiz=T)

melt.d0 <- melt(d0)
melt.d0$Sample <- factor(melt.d0$Sample,rev(c("LS174T",paste("Clone",1:14))))
  
ggplot(melt.d0,aes(y=Sample,value)) + 
  geom_col(aes(fill=variable),position=position_stack(reverse=T)) +
  geom_text(size = 3, position = position_stack(vjust = 0.5),aes(label=round(value,1)))

pG1 <- mean(d0$G1[-c(1,10)]) 
pG2 <- mean(d0$G2[-c(1,10)])

(pG1 - d0$G1) / pG1
a <- (d0$G2 - pG2) / (pG1 - pG2) 
plot(a)
round(a,2)

