#!/bin/bash
refdir=/home/mari/RSV/post-imp/1000g

cd $refdir

#Can't filter off of study data as that is in CHR:POS:BP format
#Need to convert to this format
#However multiallic snps make this v difficult
#Clean ref data then update snp annotation
#Once both datasets are in the same format, filter to use only snps that are in study data (make exlusion list from study)
#Then finally! Merge 

#Bim File Fields
#$1 = CHR, $2 = SNPID, $3 = GC, $4 = Pos $5 = Ref, $6 = Alt

#Clean dataset remove '.' snps
plink2 --bfile all_phase3 --exclude-snps . --snps-only just-acgt --make-bed --out all_phase3_qc

#79216027 varients remaining 

#remove multiallelic snps ------------------------------------
#awk reads the .bim file and outputs same-position variant IDs, then writes a list of these snps to exclude

awk '{print $2, $4}' all_phase3_qc.bim | cut -f -2 | uniq -d | cut -f -1 > multivar
awk '{print $1}' multivar > multivar_snps

#Do the same for duplicate name snps
awk '{print $2}' all_phase3_qc.bim | uniq -d >> multivar_snps
#not many dup snps detected try plink version also 

plink2 --bfile all_phase3_qc --rm-dup --make-bed --out all_phase3_qc2

cat all_phase3_qc2.rmdup.mismatch multivar_snps | uniq > exclude_list

plink2 --bfile all_phase3_qc --exclude exclude_list --make-bed --out all_phase3_qc2

#poss just try sort -k3n myFile.bim |  uniq  -f2 -D | cut -f2 > dupeSNP.txt

#Update to CHR:POS:BP -------------------------------------------

#update map name
#Make update name file
awk '{if ($1 != 0) print $2,"chr"$1":"$4":"$5":"$6}' all_phase3_qc2.bim > updateFormat.txt

echo make update name file DONE

plink2 --bfile all_phase3_qc2 --update-name updateFormat.txt --make-bed --out all_phase3_chr

head all_phase3_chr.bim

#Make list of snps from study data to keep ----------------------
#--extract

awk '{print $2}' qc_ancestry/RSV_imp_QC_remap.bim > keep_list
#only 4353 varients??

plink2 --bfile all_phase3_chr --extract keep_list --make-bed --out all_phase3_chr_toMerge

#Megrge ---------------------

plink2 --bfile ${FILEPREFIX}_QC --bmerge ${TOP}/TOP.bed ${TOP}/TOP.bim ${TOP}/TOP.fam --maf 0.1 --geno 0.1 --make-bed --out mergeTOP/mergedwTOP_test

## issue with variants at same position but different alleles - exclude these
plink2 --bfile ${FILEPREFIX}_QC --exclude mergeTOP/mergedwTOP_test-merge.missnp --make-bed --out mergeTOP/${FILEPREFIX}_forMerge

plink2 --bfile mergeTOP/${FILEPREFIX}_forMerge --bmerge ${TOP}/TOP.bed ${TOP}/TOP.bim ${TOP}/TOP.fam --maf 0.2 --geno 0.05 --make-bed --out mergeTOP/${FILEPREFIX}_mergedwTOP

echo merge DONE




