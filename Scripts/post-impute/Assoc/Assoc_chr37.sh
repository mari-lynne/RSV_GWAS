#!/bin/bash

dir=/home/mari/RSV/data/post-imp/Assoc
fuma_dir=$dir/fuma

name=RSV_Oxford_nosib
log=$dir/log

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

# Move log files

mv *.log $log

# 6) filter for ADD ------------------------------------------------------------

# First print header then filter

awk 'NR==1{print $1, $2, $4, $6, $12, $9, $10} NR!=1{if ($7 == "ADD") {print $1, $2, $4, $6, $12, $9, $10}}' \
RSV_37.RESV_total_score.glm.linear  > RSV_37_gwas.txt


#CHROM1  POS2     ID 3     REF 4    ALT 5    A1  6
#    TEST 7   OBS_CT 8  BETA 9   SE  10    T_STAT11  P 12

