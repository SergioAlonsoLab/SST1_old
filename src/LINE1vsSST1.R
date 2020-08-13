# required libraries

library(gdata)
library(ggplot2)
library(magrittr)
library(ggrepel)
library(reshape2)
library(ggsci)
library(scales)

# create figures ?

CREATE_PLOTS <- FALSE # override if plots figures should be created

# set working directory (depends on the computer) ----

try(setwd("/imppc/labs/mplab/salonso/SST1"),silent=T)
try(setwd("~/Documents/SST1/"),silent=T)
try(setwd("~/Documents/WORK/SST1/"),silent=T)

source("src/auxiliary.R") # common auxiliary functions and themes

# Graphical functions and palettes ----

# ggRegression  

ggRegression <- function(x,y,label=NULL,group=NULL,level=0.95) {
  if(is.null(label)) label <- 1:length(x)
  d0 <- data.frame(x,y)
  rownames(d0) <- label
  lm0 <- lm(y ~ x,d0)
  newX <- seq(min(x,na.rm=T),max(x,na.rm=T),l=100)
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


# short accessory functions ----






# end short accessory functions ----



# load the SST1 data from Maria/Bea ----
sst1 <- read.delim("data/SST1.MSQCR.tsv",sep="\t")

cast(sst1,Case.number ~ Type,mean)[,c(1,3,4)] %>%
  na.omit -> sst1.pairs

g1 <- ggRegression(sst1.pairs$Normal,sst1.pairs$Tumor,sst1.pairs$Case.number,level=.95) 


g1$lm %>% summary

g1$plot + xlab("SST1 RDL (log2) in Normal") + ylab("SST1 RDL (log2) in Tumor") + annotate("label",3,-2,label=expression(),parse=T)
selected <- which(sst1.pairs$Case.number %in% c(g1$high,g1$low))

ggPairedViolin(sst1.pairs[,c("Normal","Tumor")]) + sf + sc +
  geom_label_repel(data=sst1.pairs[selected,],aes(x=2,y=Tumor,label=Case.number),direction="x",nudge_x = .3) +
  xlab(NULL) + ylab("SST1 RDL (log2)") + ggtitle("SST1 RDL") 


# load LINE1 data from Sanne/Bea ----

line1 <- read.delim("data/LINE1.MSQCR.tsv",sep="\t",stringsAsFactors = F)
cast(line1,Case.number ~ Type,mean)[,c("Case.number","Normal","Tumor")] %>%
  na.omit -> line1.pairs


g2 <- ggRegression(line1.pairs$Normal,line1.pairs$Tumor,line1.pairs$Case.number)
g2$plot

ggPairedViolin(line1.pairs[,c("Normal","Tumor")]) + sc + sf + ylab("LINE1 RDL")


# common samples: samples analyzed for both SST1 and LINE1 ----

common <- merge(line1.pairs,sst1.pairs,by="Case.number",suffixes = c(".line1",".sst1"))

ggRegression(common$Tumor.line1,common$Tumor.sst1,label = common$Case.number)$plot
ggRegression(common$Normal.line1,common$Normal.sst1,label = common$Case.number)$plot
ggRegression(common$Tumor.line1,common$Tumor.sst1,label = common$Case.number)$plot

delta.line1 <- common$Tumor.line1-common$Normal.line1
delta.sst1 <- common$Tumor.sst1-common$Normal.sst1

summary(aov(delta.sst1 ~ delta.line1))
summary(lm(delta.sst1 ~ delta.line1))

ggRegression(delta.line1,delta.sst1,common$Case.number)




#----



# Analysis of nuclei size in LS174T clones ----

# Load nuclei size data

nucSize <- read.delim("data/NucleiSize.tsv",sep="\t",stringsAsFactors=F)

mSize <- tapply(nucSize$area,nucSize$sample,mean)
nucSize$sample <- factor(nucSize$sample,c("LS174T",paste("Clone",1:14)))

if(CREATE_PLOTS) png("graphs/nucleiSize.png",6,4,units = "in",res=300)
ggplot(nucSize) + geom_boxplot(aes(x=sample,y=area,fill=mSize[sample])) + 
  scale_fill_gradient(name="Mean area",low="#FFAA55",high="#CC5522") + 
  xlab(NULL) + ylab("Nuclei area (μm2)") +
  theme(axis.text.x = element_text(vjust = .5,hjust=1,angle=90)) +
  ggtitle("Nuclei area in LS174T and derived clones")
if(CREATE_PLOTS) dev.off()

aov(area ~ sample,nucSize) %>% TukeyHSD() -> tukey

nucSize$sample <- factor(nucSize$sample,c("LS174T",paste("Clone",1:14)))

table_size <- aggregate(area ~ sample,nucSize,summary)
table_size$n <- table(nucSize$sample)

write.table(table_size,file="sandbox/sizes.tsv",sep="\t",quote=F)

# Methylation of LS174T clones ----

# Load methylation data
methylation <- read.xls("data/BisulfiteClonesData.xlsx",skip=1,stringsAsFactors=F)
names(methylation)[1:2] <- c("Sample","Sequence")
coors <- methylation[1,paste("CpG",1:29,sep=".")]
methylation <- methylation[-1,]
methylation <- subset(methylation, !Sequence %in% c(".",".."))
methylation$Sample <- factor(methylation$Sample,paste("LS174T Clone",1:14)) # reorder the clones
levels(methylation$Sample) <- gsub("LS174T ","",levels(methylation$Sample))

melted0 <- melt(methylation) 

bySeq <- aggregate(value ~ Sequence + Sample,melted0,mean,na.rm=T)
bySample <- aggregate(value ~ Sample,bySeq,mean,na.rm=T)

sampleColor <- colorRamp(c("lightblue","black"))(bySample$value) %>% apply(.,1,function(x) rgb(x[1],x[2],x[3],maxColorValue = 255))

ggplot(bySeq) + geom_boxplot(aes(x=Sample,y=value),fill=sampleColor) + 
  scale_fill_meth +
  coord_cartesian(ylim=c(0,1)) + 
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=.5)) +
  xlab(NULL) + ylab("SST1 Methylation")

# Methylation vs Nuclei size ----

# There is another way to calculate average methylation, considering all sites but not by molecule


size_meth <- data.frame(Area=mSize[-1],Meth=bySample$value,Meth2=aggregate(value ~ Sample,melted0,mean)$value)

summary(lm0 <- lm(Area ~ I(1/Meth),size_meth))
summary(lm1 <- lm(Area ~ I(1/Meth2),size_meth))
summary(lm2 <- lm(Area ~ I(Meth),size_meth))
summary(lm3 <- lm(Area ~ I(Meth2),size_meth))


ggplot(size_meth,aes(x=Meth*100,y=Area)) + 
  stat_smooth(method="lm",formula=y~I(1/x)) + 
  geom_point(aes(fill=Meth),shape=21,size=4) + 
  scale_fill_meth +
  geom_label_repel(aes(label=sprintf("C%02i",1:14)),nudge_x = .015,nudge_y = runif(14,-1,1)) +
  xlab("SST1 Average Methylation (%)") +
  ylab("Average nuclei area (μm2)") +
  annotate("text",42,145,label="P=0.0035") +
  annotate("text",42,140,label="R^2==0.52",parse=T) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))



# Differences in cell cycle ----

cellCycle <- read.xls("data/CellCycleBEA.xlsx",sheet=1)
rownames(cellCycle) <- cellCycle$Sample
cellCycle <- as.matrix(cellCycle[,2:4])

colnames(cellCycle) <- c("G1","S","G2")

cellCycle <- cellCycle / rowSums(cellCycle) * 100 # scale to 100%

barplot(t(cellCycle),horiz=T,las=1)

ggplot(melt(cellCycle),aes(x=Var1,y=value,fill=Var2)) + 
  geom_col(position=position_stack(reverse=T)) + 
  coord_flip() +
  scale_fill_manual(name="Cell\ncycle\nphase",
                    labels=c("G0/1","S","G2/M"),
                    values=colorRampPalette(c("tomato3","orange"))(3)) +
  xlab(NULL) + ylab("Percentage of cells") +
  geom_text(size = 4, 
            position = position_stack(reverse = T,vjust = 0.5),
            aes(label=sprintf("%0.0f%%",value))) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


# Comparison to bisulfite sequencing in cases by Johanna ----

oldCases <- read.csv("data/SST1.csv")

melt.oldCases <- melt(oldCases)
colnames(melt.oldCases)[1] <- "Sample"
gsub("X","",melt.oldCases$Sample) -> melt.oldCases$Sample

melt.oldCases$Type <- ifelse(grepl("N",melt.oldCases$Sample),"N","T")
melt.oldCases$Patient <- gsub("[NT]","",melt.oldCases$Sample)
head(melt.oldCases)
grep(".1",melt.oldCases$Sample,fixed=T,value = T) # these are repeats

agg.oldCases <- aggregate(value ~ Type + Patient,melt.oldCases,mean) %>% cast(.,Patient ~ Type)

is.na(differ) %>% sum
dim(agg.oldCases)

# data from the excel used for NAR paper ----

allCases <- read.xls("data/SST1-all data.xls",3)

allCases$Combined.all.cases.Hypomethylation.or.No.Change

differ <- allCases$Tumor.Bisulfite.Sequencing - allCases$Normal.Bisulfite.Sequencing
sum(!is.na(differ)) # 82 samples (2 were repetitions)

bisulf.vs.msqpcr <- merge(sst1.pairs,allCases,by.x="Case.number",by.y="CASE",order=F)

names(bisulf.vs.msqpcr)

plot(I(Tumor-Normal) ~ Difference,bisulf.vs.msqpcr)
plot(I(Tumor^2) ~ Tumor.Bisulfite.Sequencing,bisulf.vs.msqpcr)
plot(Normal ~ Normal.Bisulfite.Sequencing,bisulf.vs.msqpcr)

ggplot(bisulf.vs.msqpcr) +
  aes(Tumor.Bisulfite.Sequencing,Tumor) +
  geom_smooth(method="lm") +
  geom_point(pch=21) +
  xlab("SST1 Methylation % (bisulfite sequencing)") +
  ylab("SST1 RDL (log2)") +
  geom_label_repel(aes(label=Case.number)) 

factor(bisulf.vs.msqpcr$Hypo.divided.Demethylation.and.Severe.demethylation,
       c("NC","D","D?","SD")) -> bisulf.vs.msqpcr$Classification
levels(bisulf.vs.msqpcr$Classification)[3] <- "D"

ggplot(bisulf.vs.msqpcr) +
  aes(x=Classification,
      y=Tumor-Normal) +
  geom_boxplot() +
  geom_hline(yintercept = 2,lty=2)

xtabs(~ I(Tumor - Normal >= 1.9) + Classification,bisulf.vs.msqpcr)




table(bisulf.vs.msqpcr$Hypo.divided.Demethylation.and.Severe.demethylation)

logistic1 <- glm(I(Hypo.divided.Demethylation.and.Severe.demethylation=="SD") ~ I(Tumor - Normal),bisulf.vs.msqpcr,family=binomial)
logistic2 <- glm(I(Hypo.divided.Demethylation.and.Severe.demethylation=="SD") ~ I(Tumor),bisulf.vs.msqpcr,family=binomial)

summary(logistic1)
summary(logistic2)

# prediction SD using MS-QPCR

library(ROCR)



pred1 <- prediction(bisulf.vs.msqpcr$Tumor,bisulf.vs.msqpcr$Classification=="SD")
pred2 <- prediction(bisulf.vs.msqpcr$Tumor - bisulf.vs.msqpcr$Normal,bisulf.vs.msqpcr$Classification=="SD")

plot(performance(pred1,"acc"),type="o",pch=20)
plot(performance(pred2,"acc"),type="o",pch=20)
grid()

perf1 <- performance(pred1,"tpr","fpr")
plot(perf1,type="s")

perf2 <- performance(pred2,"tpr","fpr")
plot(perf2,type="s")


x <- na.omit(data.frame(Tumor=bisulf.vs.msqpcr$Tumor.Bisulfite.Sequencing,Class=bisulf.vs.msqpcr$Hypo.divided.Demethylation.and.Severe.demethylation))
pred3 <- prediction(x$Tumor,x$Class=="SD")
plot(performance(pred3,"tpr","fpr"),type="o")
boxplot(Tumor ~ Class,x)

summary(lm(I(2^Tumor) ~ Tumor.Bisulfite.Sequencing,bisulf.vs.msqpcr))
summary(lm(I(Tumor) ~ Tumor.Bisulfite.Sequencing,bisulf.vs.msqpcr))

cor(bisulf.vs.msqpcr$Tumor,bisulf.vs.msqpcr$Tumor.Bisulfite.Sequencing,use="na")
cor.test(bisulf.vs.msqpcr$Tumor,bisulf.vs.msqpcr$Tumor.Bisulfite.Sequencing,use="na")

plot(differ,allCases$Difference)

abline(0,-1)

table(!is.na(allCases$Tumor.Bisulfite.Sequencing))

boxplot(oldCases)
summary(oldCases)


y <- melt(x) %>% dcast(.,Case.number ~ type,fun.aggregate = mean)

sst1.bisulfite <- merge(sst1.pairs,y,by="Case.number")

plot(sst1.bisulfite$T - sst1.bisulfite$N,sst1.bisulfite$Tumor - sst1.bisulfite$Normal)

# Patient's data

patients <- readRDS("data/patients.rds")

age <- patients[sst1.pairs$Case.number,"Age"]
plot(age,sst1.pairs$Tumor - sst1.pairs$Normal)

age <- patients[line1.pairs$Case.number,"Age"]
plot(age,line1.pairs$Tumor)
plot(age,line1.pairs$Normal)
plot(age,line1.pairs$Tumor - line1.pairs$Normal)

cor.test(age,line1.pairs$Tumor)
cor.test(age,line1.pairs$Normal)
cor.test(age,line1.pairs$Tumor - line1.pairs$Normal)


p53 <- patients[sst1.pairs$Case.number,"P53"]
boxplot(sst1.pairs$Tumor ~ p53)

ggplot(data.frame(sst1.pairs,p53)) + geom_boxplot(aes(x=p53,y=Tumor,fill=p53))

