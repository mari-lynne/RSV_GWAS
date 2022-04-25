#!/bin/bash'

#Split by chromosomes and output as vcf
#for i in {1..22}; do plink2 --bfile rsv_remap-updated-chr$i --chr $i --recode vcf --output-chr chrM --out rsv_remap-updated${i}; done
#echo Recode Done

#Compress

for i in {1..22}; do bcftools sort rsv_remap-updated$i.vcf -Oz -o rsv_remap-updated-chr$i.vcf.gz; done
echo RSV file compression Done

