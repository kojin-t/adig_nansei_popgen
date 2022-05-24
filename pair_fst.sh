#!bin/bash

for i in $1
do
for j in $2
do
vcftools --vcf harded5_00101.vcf \
--weir-fst-pop pops/${i}_pop \
--weir-fst-pop pops/${j}_pop \
--out pairfst_${i}_${j} > pairfst_${i}_${j}.log 2>&1
done
done

