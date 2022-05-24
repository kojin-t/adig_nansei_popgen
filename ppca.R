#!/bin/R
args <- commandArgs(trailingOnly=TRUE)
library("vcfR")
require(pcaMethods)
vcf<-read.vcfR(args)
gt<-extract.gt(vcf)
chrom <-getCHROM(vcf)
pos<-getPOS(vcf)
gt[1:10, 1:10]
gt.score<-matrix(NA, nrow(gt), ncol(gt))
gt.score[gt == "0/0" ] <- 0
gt.score[gt == "0/1"] <- 1
gt.score[gt == "1/1"] <- 2
rownames(gt.score) <- rownames(gt)
colnames(gt.score) <- colnames(gt)
gt.score <- t(gt.score)
gt.score[1:10, 1:10]

pca <- pca(gt.score, "ppca", nPcs = 10)

write.table(scores(pca), paste(args, ".ppca.result", sep=""))

write.table(summary(pca), paste(args, ".ppca.summary", sep=""))
