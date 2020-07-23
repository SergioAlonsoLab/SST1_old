# required libraries

library(gdata)
library(ggplot2)
library(magrittr)
library(ggrepel)
library(reshape2)
library(ggsci)
library(scales)


# set working directory (depends on the computer) ----

try(setwd("/imppc/labs/mplab/salonso/SST1"),silent=T)
try(setwd("~/Documents/SST1/"),silent=T)

# Graphical functions and palettes ----

# ggRegression  

ggRegression <- function(x,y,label=NULL,group=NULL,level=0.95) {
  if(is.null(label)) label <- 1:length(x)
  d0 <- data.frame(x,y)
  rownames(d0) <- label
  lm0 <- lm(y ~ x,d0)
  newX <- seq(min(x),max(x),l=100)
  pred0 <- predict.lm(lm0,newdata = data.frame(x=newX),interval = "pred",level = level)
  pred1 <- predict.lm(lm0,interval="pred",level = level,newdata = d0)
  high <- which(y > pred1[,3])
  low <- which(y < pred1[,2])
  
  gg0 <- ggplot() + aes(x=x,y=y) +
    geom_point(aes(col=group)) +
    geom_smooth(method="lm") +
    geom_line(data=data.frame(x=newX,y=pred0[,2]),lty=2,lwd=.2) +
    geom_line(data=data.frame(x=newX,y=pred0[,3]),lty=2,lwd=.2) +
    geom_label_repel(aes(x=x[high],y=y[high],label=label[high])) +
    geom_label_repel(aes(x=x[low],y=y[low],label=label[low])) 
  
  
  return(list(lm=lm0,plot=gg0,high=label[high],low=label[low]))
  
}

ggPairedViolin <- function(d0,label=NULL) {
  d1 <- melt(d0,measure.vars = 1:2,variable.name = "Group")
  Diff <- d0[,2] - d0[,1]
  g1 <- ggplot(d1,aes(x=Group,y=value)) + geom_violin(aes(fill=Group),alpha=0.1) 
  g2 <- geom_segment(aes(x=1,y=d0[,1],xend=2,yend=d0[,2],color=Diff),data=d0,lwd=0.5,inherit.aes = F) 
  g3 <- geom_point(shape=21,size=3,aes(fill=Group)) 

  return(g1+g2+g3)
  
}

# test violin
ggPairedViolin(data.frame(a=rnorm(50),b=rnorm(50))) + scale_color_gradient2(high="darkblue",low="darkred",mid = "grey")


# palettes
tnColors <- c("lightgreen","red")
sf <- scale_fill_manual(values = tnColors) 
sc <- scale_color_gradient2(mid="grey")


# end graphical functions ----


# short accesory functions ----

# find outlayers in normally distributed data

outlayers <- function(x) {
  upper <- pnorm(x,mean(x),sd(x),lower.tail = F)
  lower <- pnorm(x,mean(x),sd(x),lower.tail = T)
  data.frame(lower,upper)
}

# end short accesory functions ----



# load the SST1 data from Maria/Bea ----
sst1 <- read.xls("data/SST1_MARIA_data.xlsx",1)
sst1$TYPE <- factor(sst1$TYPE %>% gsub(" ","",.))
rownames(sst1) <- sst1$Sample

sst1$Normalized.SYBR %>% as.character() %>% as.numeric() -> sst1$Normalized.SYBR

# N vs T in SST1 ----

levels(sst1$TYPE) <- c("Normal","Tumor")

melt(sst1,measure.vars = "Normalized.SYBR",id.vars = c("Case.number","TYPE")) %>%
  dcast(., Case.number ~ TYPE,mean) %>% na.omit -> sst1.pairs
sst1.pairs[,2:3] <- log2(sst1.pairs[,2:3]) # transform to logarithmic scale

g1 <- ggRegression(sst1.pairs$Normal,sst1.pairs$Tumor,sst1.pairs$Case.number,level=.95) 


g1$lm %>% summary

g1$plot + xlab("SST1 RDL (log2) in Normal") + ylab("SST1 RDL (log2) in Tumor") + annotate("label",3,-2,label=expression(),parse=T)


selected <- which(sst1.pairs$Case.number %in% c(g1$high,g1$low))

ggPairedViolin(sst1.pairs[,c("Normal","Tumor")]) + sf + sc +
  geom_label_repel(data=sst1.pairs[selected,],aes(x=2,y=Tumor,label=Case.number),direction="x",nudge_x = .3) +
  xlab(NULL) + ylab("SST1 RDL (log2)") + ggtitle("SST1 RDL") 


# load the LINE1 data from Sanne/Bea
line1 <- read.xls("data/LINE1_summary_Bea.xlsx",2)
line1$TYPE <- factor(line1$TYPE,c("N","T"),c("Normal","Tumor"))
line1 <- subset(line1,!is.na(line1$average.RDR) & !is.na(TYPE))

melt(line1,measure.vars = "average.RDR",id.vars=c("Case.number","TYPE")) %>% 
  dcast(.,Case.number ~ TYPE,mean) %>% na.omit -> line1.pairs

line1.pairs[,2:3] <- log2(line1.pairs[,2:3])



g2 <- ggRegression(line1.pairs$Normal,line1.pairs$Tumor,line1.pairs$Case.number)
g2$plot

ggPairedViolin(line1.pairs[,c("Normal","Tumor")]) + sc + sf


# common samples: samples analyzed for both SST1 and LINE1

common <- merge(line1.pairs,sst1.pairs,by="Case.number",suffixes = c(".line1",".sst1"))

ggRegression(common$Tumor.line1,common$Tumor.sst1,label = common$Case.number)$plot
ggRegression(common$Normal.line1,common$Normal.sst1,label = common$Case.number)$plot
ggRegression(common$Tumor.line1,common$Tumor.sst1,label = common$Case.number)$plot

delta.line1 <- common$Tumor.line1-common$Normal.line1
delta.sst1 <- common$Tumor.sst1-common$Normal.sst1

ggRegression(delta.line1,delta.sst1,common$Case.number)


#----



# Analysis of nuclei size in LS174T clones ----

# Load nuclei size data

nucSize <- read.xls("data/NucleoSizeImageJ.xlsx",sheet = 3)
nucSize <- subset(nucSize,sample!="LS174T")
nucSize$sample <- factor(nucSize$sample)
levels(nucSize$sample) <- gsub("LS174_23","LS174T",levels(nucSize$sample))
levels(nucSize$sample) %>% gsub("CLONE_","Clone ",.) -> levels(nucSize$sample)
nucSize$sample <- factor(nucSize$sample,c("LS174T",paste("Clone",1:14)))
nucSize$area <- nucSize$area / 155^2 * 25^2 # 155px = 25μm 


mSize <- tapply(nucSize$area,nucSize$sample,mean)
ggplot(nucSize) + geom_boxplot(aes(x=sample,y=area,fill=mSize[sample])) + 
  scale_fill_gradient(name="Mean area",low="#FFAA55",high="#CC5522") + 
  xlab(NULL) + ylab("Area (μm2)") +
  theme(axis.text.x = element_text(angle=90)) +
  ggtitle("Nuclei area")

# Methylation of LS174T clones ----

# Load methylation data
methylation <- read.xls("data/BisulfiteClonesData.xlsx",skip=2)
methylation <- subset(methylation, Clone != "")
methylation$Clone <- factor(methylation$Clone,paste("LS174T Clone",1:14)) # reorder the clones
levels(methylation$Clone) <- gsub("LS174T ","",levels(methylation$Clone))

scale_methylation <- scale_fill_gradient(low="lightblue",high="darkblue",name="Average\nMethylation")

mMeth <- tapply(methylation$Average,methylation$Clone,mean,na.rm=T)
ggplot(methylation) + geom_boxplot(aes(x=Clone,y=Average,fill=mMeth[Clone])) + 
  scale_methylation +
  coord_cartesian(ylim=c(0,1)) + 
  theme(axis.text.x = element_text(angle=90)) +
  xlab(NULL) + ylab("Methylation")

# Methylation vs Nuclei size ----

size_meth <- data.frame(Area=mSize[-1],Meth=mMeth)

ggplot(size_meth,aes(x=Meth,y=Area)) + 
  stat_smooth(method="lm",formula=y~I(1/x)) + 
  geom_point(aes(fill=Meth),shape=21,size=4) + 
  scale_methylation +
  geom_label(aes(label=sprintf("C%02i",1:14)))


summary(lm0 <- lm(Area ~ I(1/Meth),size_meth))
summary(lm(1/Area ~ I(1/Meth),size_meth))
summary(lm(Area ~ poly(Meth,2),size_meth))
summary(lm(I(1/Area) ~ I(1/Meth),size_meth))
summary(lm(log(Area) ~ log(Meth),size_meth))


