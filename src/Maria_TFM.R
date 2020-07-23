library(gdata)
library(magrittr)

setwd("/imppc/labs/mplab/share/SST1/")

# color palette for LS174 and the clones
# this is copied from the FACS reports

clonePalette <- c("lightgrey","#96cdae","#f8dcb6","#95edfa","#e6c0f8")
names(clonePalette) <- c("grey","green","orange","blue","purple")
blue <- paste("Clone",c(1,3,11,14))
purple <- paste("Clone",c(2,4,10,12,13))
orange <- paste("Clone",5:8)
green <- paste("Clone",9)

# Load methylation data
methylation <- read.xls("data/BisulfiteClonesData.xlsx",skip=2)
methylation <- subset(methylation, Clone != "")
methylation$Clone <- factor(methylation$Clone,paste("LS174T Clone",1:14)) # reorder the clones
levels(methylation$Clone) <- gsub("LS174T ","",levels(methylation$Clone))

# Load nuclei size data

nucSize <- read.xls("data/NucleoSizeImageJ.xlsx",sheet = 3)
nucSize <- subset(nucSize,sample!="LS174T")
nucSize$sample <- factor(nucSize$sample)
levels(nucSize$sample) <- gsub("LS174_23","LS174T",levels(nucSize$sample))
levels(nucSize$sample) %>% gsub("CLONE_","Clone ",.) -> levels(nucSize$sample)
nucSize$sample <- factor(nucSize$sample,c("LS174T",paste("Clone",1:14)))
nucSize$area <- nucSize$area / 155^2 * 25^2 # 155px = 25μm



# Clones analyzed on both experiments

common <- intersect(levels(nucSize$sample),levels(methylation$Clone))

par(mfrow=c(1,1),mar=c(5,5,3,1))

meanMethylation <- tapply(methylation$Average,methylation$Clone,mean)[common]
meanSize <- tapply(nucSize$area,nucSize$sample,mean)[common]


# Plot sizes and methylation ----

layout(matrix(1:2,nrow=1),widths = c(1.75,1))

par(mar=c(5,3.5,3,1),mgp=c(2.3,.7,0))

levels(nucSize$sample) %in% green

colors <- rep(clonePalette["grey"],15)
names(colors) <- levels(nucSize$sample)
colors[blue] <- clonePalette["blue"]
colors[orange] <- clonePalette["orange"]
colors[green] <- clonePalette["green"]
colors[purple] <- clonePalette["purple"]

boxplot(area ~ sample, nucSize,las=2,xlab=NA,ylim=c(0,100),pch=19,cex=.5,ylab="Area (μm2)",col=colors)
title("Nuclei area")
#text(1:15,110,paste0("n=",table(nucSize$sample)),cex=.8)

plot(meanMethylation,meanSize,xlim=c(0.2,0.55),ylim=c(10,50),pch=NA,las=1,
     xlab="Average SST1 methylation",
     ylab="Average nuclei area (μm2)")
grid()
summary(lm0 <- lm(meanSize ~ I(1/meanMethylation)))

x0 <- seq(min(meanMethylation)-.01,max(meanMethylation)+.01,l=100)
y0 <- predict(lm0,newdata=data.frame(meanMethylation=x0),type="response",interval = "pred")

lines(x0,y0[,1],col="red",lwd=2)
lines(x0,y0[,2],col="red",lty=2)
lines(x0,y0[,3],col="red",lty=2)

abline(v=mean(meanMethylation))

points(meanMethylation,meanSize,cex=2.5,pch=21,bg=colors[common])
text(meanMethylation,meanSize,1:14,cex=.8)
text(.47,40,"R2=0.48\nP=0.0035")
title("Nuclei area vs \nSST1 methylation")

# ----

# mode sizes


sizedensities <- tapply(nucSize$area,nucSize$sample,density,from=0,to=70,bw=3)

plot(sizedensities$LS174T)

par(mfrow=c(5,3),mar=c(3,4,1.5,.5))
for(i in 1:15) {
  plot(sizedensities[[i]],main=names(sizedensities)[i],col="red",ylim=c(0,.1))
  lines(sizedensities$LS174T,col="blue")
}


mostFreqSize <- sapply(sizedensities,function(d) d$x[which.max(d$y)])
par(mfrow=c(1,1))
plot(meanSize[common],mostFreqSize[common],xlim=c(0,50),ylim=c(0,5))
grid()
abline(0,1)

# Verify that the average methylation per molecule has been correctly calculated

plot(rowMeans(methylation[,3:31],na.rm=T),methylation$Average,ylim=c(0,1),xlim=c(0,1))
abline(0,1)
grid()

# unlist methylation values per clone

values <- tapply(methylation$Average,methylation$Clone,c)

means <- sapply(values,mean)
sds <- sapply(values,sd)
medians <- sapply(values,median)
o1 <- order(means,decreasing=T)

barplot(means[o1],las=2,ylim=c(-.01,1),col="#EEFFEE") -> mids
segments(mids,means[o1],mids,means[o1]+sds[o1])
segments(mids-.2,means[o1]+sds[o1],mids+.2,means[o1]+sds[o1])
stripchart(values[o1],vertical=T,at=mids,
           add=T,method="jitter",pch=21,bg="lightblue")
abline(h=mean(means),lty=2)
title("Average methylation level in LST174 clones",ylab="Methylation level")

plot(means,sds)


# ANOVA of Karyotyped clones

type <- c(2,NA,2,NA,4,4,4,4,4,NA,2,NA,NA,2)
names(type) <- common

methylation$type <- type[methylation$Clone]


ssq <- function(aov1) {
  x <- summary(aov1)[[1]][,2] 
  names(x) <- rownames(summary(aov1)[[1]])
  return(x)
}

aov1 <- aov(Average ~ Clone,methylation)
aov2 <- aov(Average ~ type + Clone,methylation)

ssq(aov1) / sum(ssq(aov1))
ssq(aov2) / sum(ssq(aov2))

summary(aov1)
summary(aov2)

# including inferences

type2 <- c(2,2,2,2,4,4,4,4,4,2,2,4,2,2)
names(type2) <- common
methylation$type2 <- type2[methylation$Clone]

aov3 <- aov(Average ~ type + Clone,methylation)
summary(aov3)

summary(glm(factor(type2) ~ Average,methylation,family="binomial"))

# cellular cycle

facs <- read.xls("data/CellCycleBEA.xlsx")

a <- facs[,c(2,3,4)]
rownames(a) <- facs$Sample

# adjust at 100%

a <- t(apply(a,1,function(x) x/sum(x)*100))

barplot(t(a),horiz = T,las=1,col=c("lightblue","blue","blue4"))


clones2n <- paste("Clone",c(1,3,11,14))
clones4n <- paste("Clone",c(5,6,7,8))

m2nG1 <- subset(facs,Sample %in% clones2n)$X..G1 %>% mean
m2nG1.sem <- subset(facs,Sample %in% clones2n)$X..G1 %>% sd / sqrt(4)
rel.err <- m2nG1.sem / m2nG1
  
estimates <- facs$X..G1 / m2nG1
estimates.err <- estimates*rel.err

data.frame(estimates,estimates.err) * 100



