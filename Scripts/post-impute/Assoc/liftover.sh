#!/bin/bash

dir=/home/mari/RSV/data/post-imp/Assoc
fuma_dir=$dir/fuma

name=RSV_Oxford_nosib
log=$dir/log

# Notes: 21/09/22
# Run liftover.sh in dir of post-imp qc data
# Will have to rerun assoc analysis too with 37 to get plink output
# Could add any extra fuma formatting into script also 
# Call script from source so cd midway through sticks

# chr1 720239 720240 AX-32104085
# Cant do from plink output as I dont think liftover would give me ordered list

cd $dir

# 1. Print UCSC BED file ------------------------------------------------------------------
awk '{print "chr" $1, $4 -1, $4, $2 }' $name.bim > $fuma_dir/${name}_lift.bed

echo BED file DONE

cd $fuma_dir

# 2. Download (wget) chain files 38>37 
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O hg38ToHg19.over.chain.gz
        
echo Chain file download DONE


# 3. liftBed: liftOver the bed file into build 37 
liftOver ${name}_lift.bed hg38ToHg19.over.chain.gz ${name}_lift_out.bed ${name}_lift_out_unlifted.bed

echo Liftover DONE

# 4. Update data with liftover SNPs ------------------------------------------------------

# extract mapped variants snp IDs
awk '{print $4}' ${name}_lift_out.bed > ${name}_lift_map.snps
# extract updated positions
awk '{print $4, $3}' ${name}_lift_out.bed > ${name}_lift37coord.pos

echo Make SNP var and coord files DONE

#Update plink file with extracted snps and coordinates
plink2 --bfile $dir/${name} --extract ${name}_lift_map.snps --update-map ${name}_lift37coord.pos --sort-vars --make-pgen --out ${name}_remap37
#Error: Fixed-width .bed/.pgen output doesn't support sorting yet.  Generate a regular sorted .pgen first, and then reformat it
plink2 --pfile ${name}_remap37 --make-bed --out ${name}_remap37

echo Update Plink files DONE

# 5 Gchr37 Assoc Analysis ------------------------------------------------------

cd $fuma_dir

plink2 --bfile ${name}_remap37 \
--pheno ${dir}/pheno2.txt \
--pheno-name RESV_total_score \
--covar ${dir}/pca_covar2.txt \
--covar-name Sex, PC1, PC2, PC3, PC4, has_pre_ex_con, baseline_age_at_visit \
--covar-variance-standardize \
--glm \
--out RSV_37

echo Association testing DONE

# # Move log files

mv *.log $log

# 6) filter for ADD -----------------------------------------------

# First print header then filter for ADD, just keep relevant fuma columns

awk 'NR==1{print $1, $2, $4, $6, $12, $9, $10} NR!=1{if ($7 == "ADD") {print $1, $2, $4, $6, $12, $9, $10}}' \
RSV_37.RESV_total_score.glm.linear  > RSV_37_gwas.txt