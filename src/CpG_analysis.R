library(gdata)
library(magrittr)
setwd("/imppc/labs/mplab/salonso/SST1_Methylation_Maria/")

data0 <- read.xls("data/BisulfiteClonesData.xlsx",1,stringsAsFactors=F)
data0 <- data0[!data0$X %in% c("..","..."),]

meth.data <- as.matrix(data0[3:98,3:31])
meth.data <- apply(meth.data,2,as.numeric)
rownames(meth.data) <- data0$LS174T.BISULFITE.SEQUENCING.OF.SELECTED.CLONES[3:98]
colnames(meth.data) <- sprintf("CG%02i",1:29)
coor <- data0[2,3:31] %>% as.numeric
names(coor) <- colnames(meth.data)

plot(NA,xlim=range(coor)+c(0,100),
     ylim=c(nrow(meth.data)+1,0),
     ylab="Clones",
     xlab="Coordinates",
     yaxt="n"
     )
for(i in 1:nrow(meth.data)) {
  m <- meth.data[i,]
  col <- c("lightblue","red")[m+1]
  nas <- is.na(m)
  points(coor,rep(i,length(coor)),
         pch=ifelse(nas,4,22),
         
         bg=col,
         cex=ifelse(nas,.2,1))
}

BstUI <- which(diff(coor) == 2)+1
points(coor[BstUI]-1,rep(-1,length(BstUI)),pch=18,col="red")

x <- read.csv("/imppc/labs/mplab/salonso/SST1_Methylation_Maria/data/SST1allBisulfiteCLONES.csv",sep=";")




rmeans <- rowMeans(x[,-1],na.rm=T)
cor1 <- apply(x[,-1],2,cor,rmeans,use="com")

sort(cor1,decreasing = T)[1:15] %>% names -> markers

lm(rmeans ~ .,x[,markers]) %>% summary


lm(rmeans ~ .,x[,c("X600","X632","X571","X587")]) %>% summary

coor <- names(x[,-1]) %>% gsub("X","",.) %>% as.numeric

plot(coor,cor1,ylim=c(-.2,1),pch=19,xlab="Coordinate",ylab="Correlation")
segments(coor,0,coor,cor1)
abline(h=0)
rect(599,-.02,607,.02,col="red")
text(coor,cor1+.07*sign(cor1),coor,srt=90,cex=1)


cmeans <- colMeans(x[,-1],na.rm=T)

plot(coor,cmeans,ylim=c(0,1),las=1,pch=19)
segments(coor,0,coor,cmeans)
abline(h=0)
rect(599,-.02,607,.02,col="red")
text(coor,cmeans+.05,coor,srt=90,cex=1,adj=0)


