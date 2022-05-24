name <- commandArgs(trailingOnly=TRUE)[1]

library("PopGenome")

vcf<-paste("harded5_", name, "cleaned.vcf.gz", sep="")

vcf_handle<-.Call("VCF_open", vcf)

ind<-.Call("VCF_getSampleNames",vcf_handle)

chrom <- read.table(paste("harded5_", name, "cleaned.chrom", sep=""))$V1

d<-vector("numeric")

for(i in chrom){
vcf_data<-readVCF(vcf, numcols=50000, frompos=1, topos=300000000, tid=i, approx=FALSE)
vcf_data<-diversity.stats(vcf_data, pi=TRUE)
d <- c(d,vcf_data@Pi)
}

write.table(d, paste("harded5_", name, "cleaned.pi", sep=""))