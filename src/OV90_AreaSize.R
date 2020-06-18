library(gdata)
library(ggplot2)
library(magrittr)

ov90 <- read.xls("data/OV90nucleusSize.xlsx")
levels(ov90$Sample)[1] <- "OV90-CL"
ov90$Sample <- factor(ov90$Sample,c("OV90-CL",paste("OV90-C",1:16,sep="")))

means <- tapply(ov90$areas,ov90$Sample,mean)

aov1 <- aov(areas ~ Sample,ov90)
vsCL <- TukeyHSD(aov1)$Sample[1:16,]
col1 <- c("lightgreen",ifelse(vsCL[,4] < 0.01,"orange","lightblue"))

par(mar=c(8,5,2,1))
boxplot(areas ~ Sample,ov90,las=2,pch=".",log="y",col=col1)
abline(h=means[1],col="darkgreen",lwd=2)
abline(h=means[1]*c(.5,2),col="darkgreen",lty=2)

text(1:17,200,sprintf("n=%i",table(ov90$Sample)),srt=90)
title("Area size OV-90",ylab="Area")

# test with ggplot


gg1 <- ggplot(ov90,aes(x=Sample,y=areas,fill=Sample))

ggplot(ov90) + geom_boxplot(aes(Sample,areas,fill=Sample)) + 
  theme(axis.text.x = element_text(angle=90,vjust = .5)) + 
  scale_y_continuous(trans = "log10") + 
  geom_hline(yintercept = means[1]*c(0.5,1,2),lty=c(2,1,2))
  

ggplot(ov90,aes(x=Sample,y=areas,fill=Sample)) + geom_violin() + 
  geom_dotplot() +
  theme(axis.text.x = element_text(angle=90,vjust = .5)) + 
  scale_y_continuous(trans = "log10") + 
  stat_summary(fun=mean, geom="point",size=2) +
  geom_hline(yintercept = means[1]*c(0.5,1,2),lty=c(2,1,2)) 

 + geom_dotplot(aes(y=areas),binaxes='y')
  