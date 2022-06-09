library(admixtools)
library(tidyverse)
prefix='harded5_withmil'
my_f2_dir='outdir/'
extract_f2(prefix, my_f2_dir, auto_only=FALSE, maxmiss=1, format="eigenstrat")

f2_blocks = f2_from_precomp("outdir/", afprod = T)

g<-read.table("input.txt")
a<-qpgraph(f2_blocks, g, return_fstats = TRUE)
wprint(a$score)
