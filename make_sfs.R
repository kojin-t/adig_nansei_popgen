name<-commandArgs(trailingOnly=TRUE)[1]

a=read.table(paste(name,"_novar", sep=""))
b=read.table(paste(name,"_snps", sep=""))

snpnumber<-b*4*10^6/(b+a)

snpnumber<-as.numeric(snpnumber)

library("vcfR")
vcf<-read.vcfR(paste(name,"_allsites_snp_dp10_clean.vcf",sep=""))

gt<-extract.gt(vcf)
chrom <-getCHROM(vcf)
pos<-getPOS(vcf)

gt.score<-matrix(NA, nrow(gt), ncol(gt))

gt.score[gt == "0/0"] <- 0
gt.score[gt == "0/1"] <- 1
gt.score[gt == "1/1"] <- 2
gt.score<-ifelse(is.na(gt.score), 0, gt.score)
gt.score<-gt.score[sample(nrow(gt.score),snpnumber),]
test <- rowSums(gt.score)
daf_sfs <- table(test)


write.table(daf_sfs,paste(name, "_daf.sfs", sep=""), row.names=FALSE, col.names=FALSE,quote=FALSE, sep="\t")

a <- read.table(paste(name,"_daf.sfs",sep=""))

write.table(a$V2[-1],paste(name,"_sfs",sep=""), row.names=FALSE, col.names=FALSE,quote=FALSE, sep="\t")