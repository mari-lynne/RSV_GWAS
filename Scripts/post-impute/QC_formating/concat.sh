#!/bin/bash'


#Unzip
for f in *.zip; do unzip -P "V0R&rNX)Wmzcv5" -d  ~/RSV/post-imp $f; done
echo Unzip DONE

#bcftools concat

bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Ou |
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o RSV_imputed.vcf.gz

 echo Concat Done
 
 	
 	#Options:
#bcftools view -Ou -i 'R2>0.3' | you can add extra filter step here (already done in topmed)
#bcftools norm -Ou -m |
#bcftools norm -Ou -f GRCh38.fa |
#norm -m, --multiallelics -|+[snps|indels|both|any]
#split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; if SNPs and indels should be merged into a single record, specify any.
#-f, --fasta-ref FILE
#reference sequence. Supplying this option will turn on left-alignment and normalization, I think leave it normalised to topmed
