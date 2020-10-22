# Prepare data from xlsx to tsv tables 
# This script is not intended to be shared in the final version
# Only tsv tables will be shared

library(gdata)
library(magrittr)

subValue <- function(find,sub,where) {
  where[where==find] <- sub
  return(where)
}

try(setwd("~/Documents/WORK/SST1/"))
try(setwd("~/Documents/SST1/"))

# Prepare SST1 MSQPCR data ----

sst1 <- read.xls(("data/SST1_MARIA_data.xlsx"),1,stringsAsFactors=F)
sst1$TYPE <- 
  gsub(" ","",sst1$TYPE) %>%
  gsub("LS-174T","CL",.) %>% factor()

levels(sst1$TYPE) <- c("CL","Normal","Tumor")

sst1$Case.number <- 
  gsub("Cell Line","LS174T",sst1$Case.number) %>% 
  gsub("CLONE","LS174T Clone ",.)

sst1$Normalized.SYBR <- as.numeric(sst1$Normalized.SYBR) %>% log2

sst1 <- sst1[,c("Case.number","TYPE","Normalized.SYBR")]
colnames(sst1) <- c("Case.number","Type","SST1.RDL")
write.table(sst1,file="data/SST1.MSQCR.tsv",sep="\t",quote=F,row.names = F)

rm(sst1)

# Prepare LINE1 MSQPCR data ----

line1 <- read.xls("data/LINE1_summary_Bea.xlsx",2,stringsAsFactors=F)
line1 <- line1[,c("Case.number","TYPE","average.RDR")]
colnames(line1) <- c("Case.number","Type","LINE1.RDL")

line1 <- subset(line1,! is.na(LINE1.RDL))
line1$LINE1.RDL <- log2(line1$LINE1.RDL)

line1$Type <- factor(line1$Type)
levels(line1$Type) <- c("Met.Liver","Normal","Normal.Liver","Tumor")

write.table(line1,file="data/LINE1.MSQCR.tsv",sep="\t",quote=F,row.names = F)

rm(line1)

# Prepare Nuclei size data ----

nucSize <- read.xls("data/NucleoSizeImageJ.xlsx",sheet = 3)
nucSize <- subset(nucSize,sample!="LS174T")
nucSize$sample <- factor(nucSize$sample,c("LS174_23",paste("CLONE",1:14,sep="_")))
levels(nucSize$sample)[1] <- "LS174T"
levels(nucSize$sample) %>% gsub("CLONE_","Clone ",.) -> levels(nucSize$sample)
nucSize$area <- nucSize$area / 155^2 * 50^2 # 155px = 50μm in the 200x pictures # 

write.table(nucSize,file="data/NucleiSize.tsv",sep="\t",quote=F,row.names = F)

rm(nucSize)


# Prepare data from old cases. File used in NAR 2018 ----

allCases <- read.xls("data/SST1-all data.xls",3)
allCases$Classification <- factor(allCases$Hypo.divided.Demethylation.and.Severe.demethylation)
levels(allCases$Classification)[3] <- "NC"
allCases$Classification <- factor(allCases$Classification,c("NC","D","SD"))

allCases$Classification2 <- allCases$Classification
allCases$Classification2[allCases$Difference >= 0.1] <- "SD"
allCases$Classification2[allCases$Difference < 0.1] <- "D"
allCases$Classification2[allCases$Difference < 0.05] <- "NC"

table(allCases$Classification,allCases$Classification2)

names(allCases)[1] <- "Case.number"

write.table(allCases,file="data/allCases.tsv",sep="\t",quote=F,row.names = F)

rm(allCases)

# Prepare IMPPC database


imppc <- read.xls("data/IMPPC_DB.xlsx")
imppc <- imppc[,1:64]
write.table(imppc,file="data/IMPPC.tsv",sep="\t",quote=F,row.names = F)
rm(imppc)
