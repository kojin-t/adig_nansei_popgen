
for i in $each_individual
do

#quality cut and removal of adaptor sequence
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --minimum-length=25 --quality-cutoff=20 --cores=24 \
-o 0${i}_1_out.fastq \
-p 0${i}_2_out.fastq \
${i}_1.fastq \
${i}_2.fastq


#map to reference genome (Shinzato et al., 2021, Mol. Biol. Evol.)
bwa mem -M -t 12 -R "@RG\tID:foo\tSM:${i}\tPL:illumina\tLB:library1" genome.fa ${i}_1_out.fastq ${i}_2_out.fastq > ${i}.sam


#fill in mate coordinates
samtools fixmate -O bam ${i}.sam ${i}_fixmate.bam


#sort reads
samtools sort -@ 12 -O bam -o ${i}_sorted.bam -T ${i} ${i}_fixmate.bam


#produce .bai files
samtools index -@ 12 ${i}_sorted.bam


#remove duplicated reads
java -Xmx20g -jar picard.jar MarkDuplicates \
I=${i}_sorted.bam \
O=${i}_marked_duplicates.bam \
M=${i}_mark_dup_met.txt \
REMOVE_DUPLICATES=true


#produce .bai files
java -Xmx20g -jar /home/kojin/hdd1forkojin/bin/picard.jar BuildBamIndex \
I=${i}_marked_duplicates.bam


#create reailign targets
java -Xmx16g -jar GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R genome.fa \
-I ${i}_marked_duplicates.bam \
-o ${i}_marked_duplicates.bam.intervals


#realign reads
java -Xmx16g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R genome.fa \
-I ${i}_marked_duplicates.bam \
-targetIntervals ${i}_marked_duplicates.bam.intervals \
-o ${i}_realigned.bam


#first variant calls
java -Xmx3g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T UnifiedGenotyper \
-R genome.fa \
-I ${i}_realigned.bam \
-mbq 20 \
-glm BOTH \
-o ${i}_raw_variants.vcf \
-hets 0.015


#Separate SNPs and INDELs
java -Xmx3g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T SelectVariants \
-R genome.fa \
-V ${i}_raw_variants.vcf \
-selectType SNP \
-o ${i}_raw_SNP.vcf

java -Xmx3g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T SelectVariants \
-R genome.fa \
-V ${i}_raw_variants.vcf \
-selectType INDEL \
-o ${i}_raw_INDEL.vcf

#filter SNPs
java -Xmx5g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V ${i}_raw_SNP.vcf \
-G_filter "DP < 5 || DP > 20" \
-G_filterName "DP_5-20" \
-o ${i}_SNP_DPfiltered.vcf

java -Xmx5g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V ${i}_SNP_DPfiltered.vcf \
-o ${i}_DPfiltered_RankSumed_SNP.vcf \
--filterExpression 'MQRankSum < -12.5' --filterName 'LowMQRankSum' \
--filterExpression 'ReadPosRankSum < -8.0' --filterName 'LowReadPosRankSum'

java -Xmx5g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V ${i}_DPfiltered_RankSumed_SNP.vcf \
-o ${i}_DPfiltered_HDfiltered_SNP.vcf \
--filterExpression 'QD < 2.0' --filterName 'LowQD' \
--filterExpression 'FS > 60.0' --filterName 'HighFisherStrand' \
--filterExpression 'MQ < 40.0' --filterName 'LowRMSMappingQuality' \
--filterExpression 'SOR > 4.0' --filterName 'HighSOR' \
--missingValuesInExpressionsShouldEvaluateAsFailing

java -Xmx5g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T SelectVariants \
-R genome.fa \
-V ${i}_DPfiltered_HDfiltered_SNP.vcf \
-o ${i}_snp_filtered_pass.vcf \
--setFilteredGtToNocall \
--excludeFiltered \
--excludeNonVariants


#filter INDELs
java -Xmx5g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V ${i}_raw_INDEL.vcf \
-G_filter "DP < 5 || DP > 20" \
-G_filterName "DP_5-20" \
-o ${i}_INDEL_DPfiltered.vcf

java -Xmx5g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V ${i}_INDEL_DPfiltered.vcf \
-o ${i}_DPfiltered_RankSumed_INDEL.vcf \
--filterExpression 'ReadPosRankSum < -20.0' --filterName 'LowReadPosRankSum'

java -Xmx5g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V ${i}_DPfiltered_RankSumed_INDEL.vcf \
-o ${i}_DPfiltered_HDfiltered_INDEL.vcf \
--filterExpression 'QD < 2.0' --filterName 'LowQD' \
--filterExpression 'FS > 200.0' --filterName 'HighFisherStrand' \
--filterExpression 'SOR > 10.0' --filterName 'HighSOR' \

java -Xmx5g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T SelectVariants \
-R genome.fa \
-V ${i}_DPfiltered_HDfiltered_INDEL.vcf \
-o ${i}_indel_filtered_pass.vcf \
--setFilteredGtToNocall \
--excludeFiltered \
--excludeNonVariants


#combine all SNP vcf files
java -Xmx10g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T CombineVariants \
-R genome.fa \
--variant DRR099273_snp_filtered_pass.vcf \
--variant DRR099274_snp_filtered_pass.vcf \
--variant DRR099275_snp_filtered_pass.vcf \
・
・
・
--variant TU13_snp_filtered_pass.vcf \
--variant TU14_snp_filtered_pass.vcf \
--variant TU15_snp_filtered_pass.vcf \
-o all_passed_snp.vcf \
-genotypeMergeOptions UNIQUIFY


#combine all INDEL vcf files
java -Xmx10g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T CombineVariants \
-R genome.fa \
--variant 99273_indel_filtered_pass.vcf \
--variant 99274_indel_filtered_pass.vcf \
--variant 99275_indel_filtered_pass.vcf \
・
・
・
--variant TU13_indel_filtered_pass.vcf \
--variant TU14_indel_filtered_pass.vcf \
--variant TU15_indel_filtered_pass.vcf \
-o all_passed_indel.vcf \
-genotypeMergeOptions UNIQUIFY


#keep commonly-observed SNPs and INDELs 
java -Xmx8g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T SelectVariants \
-R genome.fa \
-V all_passed_snp.vcf \
-o all_known_snp.vcf \
--maxNOCALLnumber 150

java -Xmx8g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T SelectVariants \
-R genome.fa \
-V all_passed_indel.vcf \
-o all_known_indel.vcf \
--maxNOCALLnumber 150


#create realign targets
java -Xmx16g -jar GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R genome.fa \
-I ${i}_realigned.bam \
-o ${i}_realigned.bam.intervals

java -Xmx8g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R genome.fa \
-I ${i}_realigned.bam \
-targetIntervals ${i}_realigned.bam.intervals \
-known all_known_snp.vcf \
-known all_known_indel.vcf \
-o ${i}_realigned2.bam


#call variants and generate gvcf files
java -Xmx10g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R genome.fa \
-I ${i}_realigned2.bam \
-ERC GVCF \
-hets 0.015 \
-pcrModel NONE \
-nct 5 \
-o ${i}_raw_gvcf.g.vcf


#change ID
sed "/^#CHROM/{s/\t[^\t]*/\tDRR0${i}/9}" ${i}_raw_gvcf.g.vcf > ${i}_rewrite.g.vcf


#genotype all gvcf files
java -Xmx100g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R genome.fa \
--variant DRR099273_rewrite.g.vcf \
--variant DRR099274_rewrite.g.vcf \
--variant DRR099275_rewrite.g.vcf \
--variant DRR099276_rewrite.g.vcf \
・
・
・
--variant TU12_rewrite.g.vcf \
--variant TU13_rewrite.g.vcf \
--variant TU14_rewrite.g.vcf \
--variant TU15_rewrite.g.vcf \
-hets 0.015 \
-o all_combined_g.vcf


#select SNPs
java -Xmx20g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T SelectVariants \
-R genome.fa \
-V all_combined_g.vcf \
-selectType SNP \
-o all_comb_snp.vcf


#filter variants
java -Xmx25g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V all_comb_snp.vcf \
-o all_snp_RankSumed.vcf \
--filterExpression 'MQRankSum < -12.5' --filterName 'LowMQRankSum' \
--filterExpression 'ReadPosRankSum < -8.0' --filterName 'LowReadPosRankSum' 

java -Xmx25g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V all_snp_RankSumed.vcf \
-o all_snp_hdfiltered.vcf \
--filterExpression 'QD < 2.0' --filterName 'LowQD' \
--filterExpression 'FS > 60.0' --filterName 'HighFisherStrand' \
--filterExpression 'MQ < 40.0' --filterName 'LowRMSMappingQuality' \
--filterExpression 'SOR > 4.0' --filterName 'HighSOR' \
--missingValuesInExpressionsShouldEvaluateAsFailing

java -Xmx25g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R genome.fa \
-V all_snp_hdfiltered.vcf \
-G_filter "DP < 5 || DP > 50" \
-G_filterName "DP_5-50" \
-o all_snp_dpfiltered5.vcf

java -Xmx25g -jar /home/kinjyo/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
-T SelectVariants \
-R genome.fa \
-V all_snp_dpfiltered5.vcf \
--setFilteredGtToNocall \
--keepOriginalAC \
--excludeFiltered \
--maxNOCALLnumber 150 \
-o all_snp_passed5.vcf


#change FORMAT ID
sed "/^##FORMAT=<ID=AD/{s/,[^,]*/,Number=2/1}" all_selected5_snp.vcf > all_converted5_snp.vcf


#convert vcf file to plink format
bcftools annotate -Ou -x INFO all_converted5_snp.vcf |
bcftools norm -Ou -m -any |
  bcftools norm -Ou -f genome.fa |
  bcftools annotate -Ob -x ID \
    -I +'%CHROM:%POS:%REF:%ALT' \
    --vcf-idspace-to _ \
    --memory 120000 \
    --const-fid \
    --allow-extra-chr \
    --make-bed \
    --out harded5_snp

/lustre7/home/kojin-ecol/plink --memory 120000 -bfile harded5_snp --allow-extra-chr --recode --tab --out harded5_snp


#make vcf file for ppca
plink --bfile harded5_snp --memory 10000 --allow-extra-chr --recode vcf-iid --out harded5_snp


#calculate genotype missing rate per individual
/lustre7/home/kojin-ecol/plink --bfile harded5_snp --memory 40000 --missing --allow-extra-chr --out harded5_snp


#remove low-quality or misidentified samples based on ppca and missing genotype rate
/lustre7/home/kojin-ecol/plink --bfile harded5_snp --memory 10000 --allow-extra-chr --remove harded5_remove_ind.txt --make-bed --out harded5_remv


#dataset pruniing for linkage disequilibruim
/lustre7/home/kojin-ecol/plink --memory 10000 -bfile harded5_remv --allow-extra-chr --recode --out harded5_rem
/lustre7/home/kojin-ecol/plink --file harded5_remv --memory 10000 --allow-extra-chr --indep-pairwise 10 5 0.2 --out range_comp
/lustre7/home/kojin-ecol/plink --bfile harded5_remv --memory 10000 --extract range_comp.prune.in --genome --allow-extra-chr --make-bed --out harded5_prune


#dataset quality controlling for mapping allele frequency, missing genotype rate, and HWE
/lustre7/home/kojin-ecol/plink --bfile harded5_prune --memory 10000 --maf 0.01 --hwe 0.001 --geno 0.1 --allow-extra-chr --make-bed --out harded5_clean

#generate vcf file from plink file
/lustre7/home/kojin-ecol/plink --bfile harded5_clean --memory 10000 --allow-extra-chr --recode vcf-iid --out harded5_clean
