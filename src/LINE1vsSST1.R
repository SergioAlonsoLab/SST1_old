# required libraries

library(gdata)
library(ggplot2)
library(magrittr)
library(ggrepel)

# load the LINE1 data from Sanne/Bea
line1 <- read.xls("data/LINE1_summary_Bea.xlsx",2)
line1 <- subset(line1,!is.na(line1$average.RDR))
line1$Sample %>% .[duplicated(.)] -> dups
dups
# to simplify, duplicated analyses were removed

line1 <- line1[-which(duplicated(line1$Sample)),]
rownames(line1) <- line1$Sample

# load the SST1 data from Maria/Bea
sst1 <- read.xls("data/SST1_MARIA_data.xlsx",1)
sst1$TYPE <- factor(sst1$TYPE %>% gsub(" ","",.))
rownames(sst1) <- sst1$Sample


# common samples

samples <- intersect(sst1$Sample,line1$Sample)

# create a dataframe with the sst1 and line1 values

data0 <- sst1[samples,]
data0$line1 <- line1[samples,"average.RDR"]
data0$line1.error <- line1[samples,"Error.propagation.1"]

normals <- subset(data0,TYPE=="N")
tumors <- subset(data0,TYPE=="T")

all(normals$Case.number == tumors$Case.number) # check they are in the same order

# ggRegression

ggRegression <- function(x,y,label=NULL,group=NULL) {
  if(is.null(label)) label <- 1:length(x)
  d0 <- data.frame(x,y)
  rownames(d0) <- label
  lm0 <- lm(y ~ x,d0)
  newX <- seq(min(x),max(x),l=100)
  pred0 <- predict.lm(lm0,newdata = data.frame(x=newX),interval = "pred")
  pred1 <- predict.lm(lm0,interval="pred")
  high <- which(y > pred1[,3])
  low <- which(y < pred1[,2])
  
  gg0 <- ggplot() + aes(x=x,y=y) +
    geom_point(aes(col=group)) +
    geom_smooth(method="lm") +
    geom_line(data=data.frame(x=newX,y=pred0[,2]),lty=2,lwd=.2) +
    geom_line(data=data.frame(x=newX,y=pred0[,3]),lty=2,lwd=.2) +
    geom_label_repel(aes(x=x[high],y=y[high],label=label[high])) +
    geom_label_repel(aes(x=x[low],y=y[low],label=label[low])) 
    
  
  return(list(lm=lm0,plot=gg0))
  
}

g1 <- ggRegression(log2(data0$Normalized.SYBR),log2(data0$Normalized.HEX),label=data0$Sample,group = data0$TYPE) 
g1$plot + xlab("SST1 (SYBR)") + ylab("SST1 (HEX)")

g2 <- ggRegression(log2(normals$Normalized.SYBR),log2(tumors$Normalized.SYBR),label=tumors$Sample)
g2$plot + xlab("SST1 (SYBR) in Normal") + ylab("SST1 (SYBR) in Tumor")

g3 <- ggRegression(log2(normals$Normalized.HEX),log2(tumors$Normalized.HEX),label=tumors$Sample)
g3$plot + xlab("SST1 (HEX) in Normal") + ylab("SST1 (HEX) in Tumor")

g4 <- ggRegression(log2(tumors$line1),log2(tumors$Normalized.SYBR),label=tumors$Sample)
g4$plot + xlab("LINE1 RDL") + ylab("SST1 (SYBR) RDL")

g5 <- ggRegression(log2(tumors$line1),log2(tumors$Normalized.HEX),label=tumors$Sample)
g5$plot + xlab("LINE1 RDL") + ylab("SST1 (HEX) RDL")

g6 <- ggRegression(log2(data0$line1),log2(data0$Normalized.SYBR),label=data0$Sample,group=data0$TYPE)
g6$plot + xlab("LINE1 RDL") + ylab("SST1 (SYBR) RDL")

g7 <- ggRegression(log2(data0$line1),log2(data0$Normalized.HEX),label=data0$Sample,group=data0$TYPE)
g7$plot + xlab("LINE1 RDL") + ylab("SST1 (HEX) RDL")

# differece between HEX and SYBR

lm0 <- lm(log2(Normalized.HEX) ~ log2(Normalized.SYBR),data0)
summary(lm0)
newX <- seq(0.001,1,l=100)
pred0 <- predict(lm0,newdata = data.frame(Normalized.SYBR=newX),interval = "pred")

ggplot(data0) + 
  aes(x=log2(Normalized.SYBR),y=log2(Normalized.HEX)) +
  geom_point(aes(col=TYPE)) +
  geom_smooth(method="lm") +
  geom_line(data=data.frame(),aes(x=log2(newX),y=pred0[,3]))



# create a general ggplot of data0

gg1 <- ggplot(data0)

gg1 + aes(x=log2(Normalized.HEX),y=log2(Normalized.SYBR)) + geom_point() + geom_smooth(method="lm")
gg1 + geom_point(aes(log2(Normalized.HEX),log2(line1))) 
gg1 + aes(x=log2(Normalized.HEX),y=log2(line1)) + geom_point() + geom_smooth(method="lm")

# boxplot of line1 and sst1 in normal and tumors

gg1 + aes(x=TYPE,y=log2(Normalized.HEX)) + geom_violin(aes(fill=TYPE))
gg1 + aes(x=TYPE,y=log2(line1)) + geom_violin(aes(fill=TYPE))

tumors <- subset(data0,TYPE=="T")
normals <- subset(data0,TYPE=="N")

rownames(tumors) %>% gsub("T","N",.) == rownames(normals) # checking they are paired

t.test(log2(tumors$Normalized.SYBR),log2(normals$Normalized.SYBR),paired=T) # no difference 
t.test(log2(tumors$line1),log2(normals$line1),paired=T)

# difference tumors vs normals

plot(log2(tumors$line1)-log2(normals$line1),log2(tumors$Normalized.HEX)-log2(normals$Normalized.HEX))

gg1 + geom_density(aes(Normalized.HEX,fill=TYPE)) + scale_x_log10()
gg1 + geom_density(aes(line1,fill=TYPE)) + scale_x_log10()
gg1 + geom_point(aes(y=Normalized.HEX,x=line1,col=TYPE)) + scale_x_log10() + scale_y_log10()

barplot(data0$Normalized.SYBR,col=ifelse(data0$TYPE=="T","red","lightblue"))

ggplot() +
  aes(x=normals$Normalized.HEX,y=tumors$Normalized.HEX) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline()


# SST1 in Normals vs Tumors

sst1_N_T <- data.frame(Patient = tumors$Case.number,
                       Normal = log2(normals$Normalized.HEX),
                       Tumor = log2(tumors$Normalized.HEX))


ggplot(data0) + geom_violin(aes(x=TYPE,y=log2(Normalized.HEX))) +
  geom_point(aes(x=TYPE,y=log2(Normalized.HEX),col=TYPE),size=3) 


selected <- which(abs(log2(normals$Normalized.HEX)-log2(tumors$Normalized.HEX)) > 1.5 )
ggplot() + 
  geom_violin(aes(x="N",y=normals$Normalized.HEX)) +
  geom_violin(aes(x="T",y=tumors$Normalized.HEX)) + 
  geom_segment(aes(x="N",y=normals$Normalized.HEX,xend="T",yend=tumors$Normalized.HEX),lwd=.1) +
  geom_point(aes(x=data0$TYPE,y=data0$Normalized.HEX,col=data0$TYPE),size=3) +
  geom_label_repel(aes(x="T",y=tumors$Normalized.HEX[selected],label=tumors$Sample[selected]),
                   nudge_x=15,segment.size=.1) +
  geom_label_repel(aes(x="N",y=normals$Normalized.HEX[selected],label=normals$Sample[selected]),
                   nudge_x=-15,segment.size=.1) +
  scale_y_log10() +
  xlab("Sample Type") +
  ylab("SST1 RDL (HEX)") +
  ggtitle("SST1 RDL by Sample Type") +
  theme(legend.position = "none",plot.title = element_text(hjust=0.5)) 
  

selected <- which(abs(log2(normals$Normalized.SYBR)-log2(tumors$Normalized.SYBR)) > 1.5 )
ggplot() + 
  geom_violin(aes(x="N",y=normals$Normalized.SYBR)) +
  geom_violin(aes(x="T",y=tumors$Normalized.SYBR)) + 
  geom_segment(aes(x="N",y=normals$Normalized.SYBR,xend="T",yend=tumors$Normalized.SYBR),lwd=.1) +
  geom_point(aes(x=data0$TYPE,y=data0$Normalized.SYBR,col=data0$TYPE),size=3) +
  geom_label_repel(aes(x="T",y=tumors$Normalized.SYBR[selected],label=tumors$Sample[selected]),
                   nudge_x=15,segment.size=.1) +
  geom_label_repel(aes(x="N",y=normals$Normalized.SYBR[selected],label=normals$Sample[selected]),
                   nudge_x=-15,segment.size=.1) +
  scale_y_log10() +
  xlab("Sample Type") +
  ylab("SST1 RDL (SYBR)") +
  ggtitle("SST1 RDL by Sample Type") +
  theme(legend.position = "none",plot.title = element_text(hjust=0.5)) 




plotNvsT <- function(foo,title="") {
  
  foo$log2ratio <- log2(foo$y/foo$x)
  foo$class <- cut(foo$log2ratio,c(-Inf,-1.5,1.5,Inf))
  levels(foo$class) <- c("Hyper","NC","Hypo")
  selected <- abs(foo$log2ratio) > 1.5 
  ggplot(data = foo) + 
    geom_point(aes(log2(x),log2(y),col=class),size=2) +
    geom_abline(intercept = c(1.5,0,-1.5), lty=c(2,1,2)) +
    geom_label_repel(data=foo[selected,],aes(log2(x),log2(y),label=caseNumber)) +
    xlab("RDL in Normal Tissue (log2)") +
    ylab("RDL in Tumor Tissue (log2)") +
    ggtitle(title)
  
}


foo <- data.frame(x=normals$Normalized.SYBR,y=tumors$Normalized.SYBR,caseNumber=tumors$Case.number)
plotNvsT(foo)

foo <- data.frame(x=normals$Normalized.HEX,y=tumors$Normalized.HEX,caseNumber=tumors$Case.number)
plotNvsT(foo)


ggplot(data0) +
  geom_point(aes(x=Normalized.HEX,y=Normalized.SYBR,col=TYPE)) +
  scale_x_log10() + scale_y_log10()




ggplot(data0) + geom_dotplot(aes(x=TYPE,y=log2(Normalized.HEX)),binwidth = 1/5)
ggplot(data0) + aes(x=TYPE,y=log2(Normalized.HEX)) + geom_violin() + geom_boxplot(width=.2)


# scrapbook

gg1 <- ggplot(subset(line1,TYPE %in% c("N","T")),aes(TYPE,average.RDR,Error.propagation.1))

gg1 + geom_boxplot(aes(TYPE,log10(average.RDR)))

myscale <- c("darkgreen","darkred")
vp1 <- gg1 + aes(TYPE,log2(average.RDR),fill=TYPE) + geom_violin()

vp1 + scale_fill_brewer(palette=7) + labs(y="log2 RDL",x="Sample Type") + theme_light()



gg2 <- ggplot(sst1)

ggplot(data.frame(x=1:10,y=10:1)) +  geom_density(aes(x))

ggplot(sst1) + aes(log2(sst1$Normalized.SYBR),log2(sst1$Normalized.HEX)) + geom_point() + geom_smooth(method="lm")

shapiro.test(sst1$Normalized.SYBR)
shapiro.test(log2(sst1$Normalized.SYBR))

ggplot(sst1) + geom_density(aes(log2(sst1$Normalized.SYBR)),fill="lightblue",alpha=.6) +
  geom_density(aes(log2(sst1$Normalized.HEX)),fill="pink",alpha=.6)
