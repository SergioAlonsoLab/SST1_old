# THIS SCRIPT IS FOR TESTING ONLY
# NOT TO BE INCLUDED IN THE FINAL VERSION OF THE MANUSCRIPT

library(gdata)
library(ggplot2)
library(magrittr)
library(igraph)
library(scales)
library(mixtools)
library(reshape)

try(setwd("~/Documents/WORK/SST1/"),silent=T)


rotate.x.labels <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

nucleiSize <- read.xls("data/NucleoSizeImageJ.xlsx",3)
nucleiSize <- subset(nucleiSize,sample!="LS174T")
nucleiSize$sample <- factor(nucleiSize$sample,c("LS174_23",paste("CLONE",1:14,sep="_")))
levels(nucleiSize$sample)[1] <- "LS174T"

nucleiSize$area <- nucleiSize$area * (25 / 155)^2

tapply(nucleiSize$area,nucleiSize$sample,shapiro.test)

png("graphs/nucleiSize.png",600,600)
ggplot(nucleiSize,aes(x=sample,y=area)) + geom_boxplot(fill="lightblue") + 
  rotate.x.labels + xlab(NULL) + ylab("Nuclei area (µm2)") + 
  theme(text=element_text(size = 30))
dev.off()

aov(area ~ sample,nucleiSize) %>% TukeyHSD() -> tukey

tukey

strsplit(rownames(tukey$sample),split="-") %>% 
  unlist %>% 
  matrix(.,byrow=T,ncol=2) -> edgelist

g0 <- graph_from_edgelist(edgelist,directed = F)
E(g0)$weight <- 1 / abs(tukey$sample[,1])
E(g0)$curved <- .25
V(g0)$name <- gsub("CLONE_","",V(g0)$name)
V(g0)$size <- 20
V(g0)$color <- "lightblue"
set.seed(123)
plot(g0,layout=layout_with_fr(g0))

meanSize <- tapply(nucleiSize$area,nucleiSize$sample,mean) 
meanSize / meanSize[1] * 100

areaBySample <- tapply(nucleiSize$area,nucleiSize$sample,c)

sapply(areaBySample,summary) %>% t %>% round(.,1) %>% write.table(.,file="sandbox/sizes.txt",sep="\t",quote=F)

mm <- normalmixEM(areaBySample$LS174T) 

ggplot(nucleiSize,aes(x=area)) + geom_density(fill="red",alpha=.2) + facet_wrap(sample ~ .)
