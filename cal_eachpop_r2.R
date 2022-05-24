name <- commandArgs(trailingOnly=TRUE)[1]

library(windowscanr)

a<-read.table(paste(name, "_ld_table.txt", sep=""), header = T)

b<-data.frame(a$CHR_B-a$CHR_A, a$SNP_B)

c<-b[order(b$a.CHR_B...a.CHR_A),]

colnames(c)<-c("base", "r2")

d<-pos_win<-winScan(x=c,
                 position="base",
                 values="r2",
                 win_size=100,
                 win_step=100)

pdf(paste(name, "_5_snp_cumr2.pdf", sep=""), height=10, width=20)

plot(pos_win$win_mid,pos_win$r2_mean,type = "l",xlim=c(0,100000), ylim=c(0.4,0.8))

dev.off()
