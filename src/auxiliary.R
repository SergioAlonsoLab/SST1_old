# ggplot themes ----

no.x <- theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
# methylation scale ----
scale_fill_meth <- scale_fill_gradient(low="lightblue",
                                        high="black",
                                        na.value = "white",
                                        name="\u03b2-value") 

# find outlayers in normally distributed data ----
outlayers <- function(x) {
  upper <- pnorm(x,mean(x),sd(x),lower.tail = F)
  lower <- pnorm(x,mean(x),sd(x),lower.tail = T)
  data.frame(lower,upper)
}

barCode <- function(cancertype="CRC",patient,sampletype=0,sample=1) {
  sprintf("%s-%06d-%02d-%02d",cancertype,patient,sampletype,sample)
}