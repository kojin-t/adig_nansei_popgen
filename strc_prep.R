args <- commandArgs(trailingOnly=TRUE)
library("vcfR")
library(RColorBrewer)
require(pcaMethods)
vcf<-read.vcfR(args)
gt<-extract.gt(vcf)
chrom <-getCHROM(vcf)
pos<-getPOS(vcf)
"gt"
gt[1:10, 1:10]
gt.score<-matrix(NA, nrow(gt), ncol(gt))

gt.score[gt == "0/0"] <- 0
gt.score[gt == "0/1"] <- 1
gt.score[gt == "1/1"] <- 2

gt.score<-ifelse(is.na(gt.score),-9,gt.score)
gt.score<-t(gt.score)

library(dplyr)
gt.contra<-gt.score[,sample(ncol(gt.score),150000)]

write.table(gt.contra, "gt.contra")