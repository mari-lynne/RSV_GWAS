#!/bin/bash

#Remove reference duplicated snps -------------------------------------------------
#awk script reads the Ref data .bim file and outputs same-position variant IDs
#Then writes a list of these snps to exclude using plink

#Remove multi-allelic snps ($4=Pos)
awk '{print $2, $4}' $qcdir/$refname.no_acgt.bim | cut -f -2 | uniq -d | cut -f -1 > multivar
awk '{print $1}' multivar > multivar_snps 

#Remove duplicate name snps ($2=SNPID)
awk '{print $2}' $qcdir/$refname.no_acgt.bim | uniq -d >> multivar_snps #add to list

#Try Plink Rmdup function also
plink2 --bfile $qcdir/$refname.no_acgt --rm-dup --out $qcdir/$refname.no_acgt2

#Combine (Cat) bash multivar snps and plink duplicate snps, save as exclude list
cat $qcdir/$refname.no_acgt2.rmdup.mismatch multivar_snps | uniq > exclude_list

#Exclude from ref dataset
plink2 --bfile $qcdir/$refname.no_acgt --exclude exclude_list --make-bed --out $qcdir/$refname.rm_dup


#removing by position will exclude too many as also need to combine with chromosome data from CHR1

#Remove multi-allelic snps ($4=Pos)

#Remove multi-allelic snps ($4=Pos)
awk '{print $2, $4}' $qcdir/$refname.no_acgt.bim | cut -f -2 | uniq -d | cut -f -1 > multivar
awk '{print $1}' multivar > multivar_snps

#change to this 

awk '{print $1, $2, $4}' $qcdir/$refname.no_acgt.bim | cut -f 1,3 | uniq -d | > multivar
awk '{print $1}' multivar > multivar_snps 


#also test exclude-all option
