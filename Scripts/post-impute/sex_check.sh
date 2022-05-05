#!/bin/bash

#Sex check redo
#Previous run still has 65 females with very high F scores (homozygous), therefore being identified as males

#Double check variant quality (should have been qc'd)
#Redo steps
#Remove XY region - Just select XY chrs
#Add QC thresholds
#LD prune
#Then run split-x 37
#Then remove (set missing hh)
#Then run sex-check with y counts

cd /home/mari/RSV/data/plink

plink2 \
--bfile rsv_QC2 \
--keep study_samples.txt \
--update-sex update_sex.txt \
--make-bed \
--out rsv_sex

plink \
--bfile rsv_sex \
--chr 23-24 \
--maf 0.01 \
--geno 0.05 \
--mind 0.05 \
--hwe 0.000001 \
--make-bed \
--out rsv_sex_chr
#6090 variants 28
#Warning: 283616 het. haploid genotypes present (see rsv_sex_chr.hh )
#Warning: Nonmissing nonmale Y chromosome genotype(s) present


plink \
--bfile rsv_sex_chr \
--exclude exclude_snps_xy.txt \
--make-bed \
--out rsv_sex_qc

#Warning: 69485 het. haploid genotypes present (see rsv_sex_qc.hh ); many
#Total genotyping rate is 0.997583.
#5804 variants and 318 people pass filters and QC.

plink \
--bfile rsv_sex_qc \
--split-x b37 \
--make-bed \
--out rsv_sex_chr

#43 chromosome codes changed, is that 43 variants?
#hh genotypes stayed the same so not an autosome region problem

plink \
--bfile rsv_sex_chr \
--set-hh-missing \
--indep-pairphase 20000 2000 0.5 \
--out rsv_sex_LD

#Warning: 66700 het. haploid genotypes present (see rsv_sex_LD.hh )
#5804 variants and 318 people pass filters and QC.
#Pruned 2065 variants from chromosome 23, leaving 3696.
#Pruned 14 variants from chromosome 25, leaving 29.
#Pruning complete.  2079 of 5804 variants removed.


plink \
--bfile rsv_sex_chr \
--extract rsv_sex_LD.prune.in \
--make-bed \
--out rsv_sex_clean
#Warning: 42368 het. haploid genotypes present (see rsv_sex_clean.hh ); many
#3725 variants and 318 people pass filters and QC.

plink \
--bfile rsv_sex_clean \
--set-hh-missing \
--make-bed \
--out rsv_sex_hh

#65 mismatches with het-haplotypes removed
#145 mismatches without

plink \
--bfile rsv_sex_hh \
--check-sex ycount  \
--out rsv_redo

#Still 65 errors

#Errors
#"<number> het. haploid genotypes present (see plink.hh )."
#This is usually caused by male heterozygous calls in the X chromosome pseudo-autosomal region. Check the variants named in the .hh file; if they are all near the beginning or end of the X chromosome, --split-x should solve the problem.
