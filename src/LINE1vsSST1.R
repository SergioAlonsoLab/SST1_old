# required libraries

library(gdata)
library(ggplot2)
library(magrittr)

# load the LINE1 data from Sanne/Bea
line1 <- read.xls("data/LINE1_summary_Bea.xlsx",2)

# load the SST1 data from Maria/Bea
sst1 <- read.xls("data/SST1_MARIA_data.xlsx",1)

sst1$Normalized.SYBR %>% as.character %>% as.numeric -> sst1$Normalized.SYBR
sst1$Normalized.HEX %>% as.character %>% as.numeric -> sst1$Normalized.HEX



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
