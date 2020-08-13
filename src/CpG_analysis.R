library(gdata)
library(magrittr)
library(ggplot2)
library(ggrepel)

try(setwd("/imppc/labs/mplab/salonso/SST1_Methylation_Maria/"),silent=T)
try(setwd("~/Documents/WORK/SST1/"),silent=T)

CREATE_PLOTS <- TRUE

source("src/auxiliary.R") # common auxiliary functions and themes

# ----

data0 <- read.xls("data/BisulfiteClonesData.xlsx",1,stringsAsFactors=F,skip=1)
names(data0)[1] <- "Sample"

cpgs <- paste("CpG",1:29,sep=".")
coors <- unlist(data0[1,cpgs])

data0 <- data0[!data0$Label %in% c("..","..."),] # remove spacers

melted0 <- melt(data0[-1,])
melted0$Sample <- factor(melted0$Sample,paste("LS174T Clone",1:14))

head(melted0)

bySample <- aggregate(value ~ Sample + variable,melted0,mean,na.rm=T)
names(bySample) <- c("Sample","CpG","value")
bySample$Sample <- factor(bySample$Sample,paste("LS174T Clone",1:14))

if(CREATE_PLOTS) png("graphs/bisulfiteSeq1.png",8,4,units="in",res=300)
ggplot(bySample,aes(x=coors[CpG],y=Sample)) +
  geom_hline(aes(yintercept=Sample)) +
  geom_point(aes(fill=value),pch=21,size=4) +
  scale_fill_meth +
  xlab("CpG coordinates (bp)")
if(CREATE_PLOTS) dev.off()

if(CREATE_PLOTS) png("graphs/bisulfiteSeq2.png",8,8,units="in",res=300)
ggplot(melted0,aes(x=variable,y=Label)) +
  geom_hline(aes(yintercept=Label)) +
  geom_point(aes(fill=value),pch=21,size=2) +
  scale_fill_meth +
  facet_wrap(~ Sample,ncol = 3,scales = "free") +
  theme_bw() +
  ylab("") +
  no.x
if(CREATE_PLOTS) dev.off()

aggregate(value ~ Sample,melted0,mean)

# global average, not by molecule ----

gmeans <- aggregate(value ~ Sample,melted0,mean,na.rm=T)

# average by molecule ----

mol.means <- aggregate(value ~ Label + Sample,melted0,mean,na.rm=T)
levels(mol.means$Sample) <- gsub("LS174T ","",levels(mol.means$Sample))

ggplot(mol.means) +
  geom_boxplot(aes(Sample,value)) +
  theme(axis.text.x = element_text(angle=90,vjust=.5,hjust = 1)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_hline(yintercept=.35) 

aov(value ~ Sample,mol.means) %>% summary

aggregate(value ~ Label + Sample,melted0,mean,na.rm=T) %>% aggregate(value ~ Sample,.,summary)

# Correlation of individual CpGs with the SST1 mean ----

meth.data <- as.matrix(data0[-1,cpgs])

mmeth <- rowMeans(meth.data,na.rm=T)

plot(gmeans$value,mmeth)

plot(mmeth,data0$AVERAGE[-1]) # check that the means were correctly calculated
mcor <- apply(meth.data,2,cor,mmeth,use="na")
mcor <- data.frame(coors,mcor,color="#FFFFFFBB")
mcor$color[mcor$coors %in% c(496,506,607,609)] <- "lightblue"

if(CREATE_PLOTS) png("graphs/correlation.png",6,4,units="in",res=300)
ggplot(mcor,aes(coors,mcor)) +
  geom_segment(xend=coors,yend=0,lty=3) +
  geom_hline(yintercept = 0) +
  xlab("CpG site coordinate (bp)") +
  ylab("Correlation with average methylation") +
  scale_y_continuous(limits = c(-.1,1)) +
  geom_label(label=1:29,fill=mcor$color) +
  geom_label(x=506-10,y=0,label="SST1-UF",size=2) +
  geom_label(x=607+10,y=0,label="SST1-UR",size=2)
if(CREATE_PLOTS) dev.off()



bySample2 <- cast(bySample,formula = Sample ~ CpG) 
rownames(bySample2) <- bySample2$Sample
bySample2$Sample <- NULL

hclust(dist(bySample2)) %>% plot

el <- dist(bySample2,upper = F) %>% as.matrix %>% melt

el$Var1 <- gsub("LS174T Clone ","",el$Var1)
el$Var2 <- gsub("LS174T Clone ","",el$Var2)

g0 <- graph_from_edgelist(as.matrix(el[,1:2]),directed=F)
E(g0)$weight <- (2-el$value)^2
plot(g0,layout=layout_with_kk(g0))

# -----

g0 <- graph_from_literal(A-B,A-C,B-C)
E(g0)$weight <- c(20,1,.1)
plot(g0,layout=layout_with_kk(g0))  

