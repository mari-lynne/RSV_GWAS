#!/bin/bash

## REQUIRES the following software
# plink,  gcta

## REQUIRES the following variables in config file
# RAWDATADIR, FILEPREFIX, TOP

## INPUT
# ${FILEPREFIX}_QCd # binary plink files following prelim QC

## OUTPUT
# mergeTOP/${FILEPREFIX}_mergedwTOP # variants merged with 100 genomes and filtered to common, shared variants
# mergeTOP/${FILEPREFIX}_mergedwTOP.pca # pca for sample and 1000 genome combined

PROCESSDIR=/home/mari/RSV/post-imp/ #Run script in /home/mari/RSV/post-imp/ (parent folder of all the data)
FILEPREFIX=RSV
TOP=/home/mari/RSV/Plink/topmed #file location #TOPMED38.vcf.gz

cd ${PROCESSDIR}
mkdir -p mergeTOP

# Change variant ids to chr:bp
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' ${FILEPREFIX}_QCd.bim > updateToTOPFormat.txt
#plink2 --bfile ${FILEPREFIX}_QC --update-name updateToTOPFormat.txt --make-bed --out ${FILEPREFIX}_QC

# Convert Topmed data to plink bfile

plink2 --vcf ${TOP}TOPMED38.vcf.gz --out ${TOP}/TOP
echo TOPMED convert vcf to PLINK DONE 

# First merge with topmed data and filter variants to those in common
# Need to test initially in case of error with triallelic variants
plink2 --bfile ${FILEPREFIX}_QC --bmerge ${TOP}/TOP.bed ${TOP}/TOP.bim ${TOP}/TOP.fam --maf 0.1 --geno 0.1 --make-bed --out mergeTOP/mergedwTOP_test

## Issue with variants at same position but different alleles - exclude these
plink2 --bfile ${FILEPREFIX}_QC --exclude mergeTOP/mergedwTOP_test-merge.missnp --make-bed --out mergeTOP/${FILEPREFIX}_forMerge

plink2 --bfile mergeTOP/${FILEPREFIX}_forMerge --bmerge ${TOP}/TOP.bed ${TOP}/TOP.bim ${TOP}/TOP.fam --maf 0.2 --geno 0.05 --make-bed --out mergeTOP/${FILEPREFIX}_mergedwTOP

echo merge DONE

# LD prune
plink2 --bfile mergeTOP/${FILEPREFIX}_mergedwTOP --indep 50 5 1.5 --out mergeTOP/${FILEPREFIX}_mergedwTOP.ld
plink2 --bfile mergeTOP/${FILEPREFIX}_mergedwTOP --extract mergeTOP/${FILEPREFIX}_mergedwTOP.ld.prune.in --make-bed --out mergeTOP/${FILEPREFIX}_mergedwTOP.ld.prune

echo pruning DONE

rm mergeTOP/${FILEPREFIX}_forMerge
rm mergeTOP/mergedwTOP_test*

echo clean files DONE

# Use GCTA to calc PCs
#${GCTA}/gcta64 --bfile mergeTOP/${FILEPREFIX}_mergedwTOP.ld.prune --make-grm-bin --autosome --out mergeTOP/${FILEPREFIX}_mergedwTOP
#${GCTA}/gcta64 --grm mergeTOP/${FILEPREFIX}_mergedwTOP --pca --out mergeTOP/${FILEPREFIX}_mergedwTOP.pca

# rm mergeTOP/${FILEPREFIX}_mergedwTOP*grm*

# plot PCs
#Rscript ${SCRIPTDIR}/utilitys/plotEthnicity.r ${PROCESSDIR}/mergeTOP/ ${FILEPREFIX} ${TOP} 
