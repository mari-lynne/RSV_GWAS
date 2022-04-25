#!/bin/bash

#Run in terminal as source

refdir=/home/mari/RSV/post-imp/1000g
mkdir -p $refdir/plink_log
log=$refdir/plink_log

cd $refdir

#Download 1000g data -------------------------------------------------------

pgen=https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1
sample=https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1
wget $pgen
mv 'all_phase3.pgen.zst?dl=1' all_phase3.pgen.zst
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
wget $pvar
mv 'all_phase3.pvar.zst?dl=1' all_phase3.pvar.zst
wget $sample
mv 'phase3_corrected.psam?dl=1' all_phase3.psam

echo 1000g data download DONE

plink2 --pfile $refdir/all_phase3 vzs --max-alleles 2 --make-bed --out $refdir/all_phase3
mv $refdir/all_phase3.log $log


#Combine with QC data-----------------------------------------------------

mkdir -p $refdir/qc_ancestry

qcdir=$refdir/qc_ancestry
name=RSV_imp_QC
log=$refdir/plink_log

cp [RSV_imp_QC]* /home/mari/RSV/post-imp/qc_ancestry

cd $qcdir

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $qcdir/$name.bim  > \
    $qcdir/$name.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $refdir/$refname.bim  > \
    $qcdir/$refname.ac_gt_snps
   
plink --bfile  $refdir/$refname \
      --exclude $qcdir/$refname.ac_gt_snps \
      --make-bed \
      --out $qcdir/$refname.no_ac_gt_snps


plink --bfile  $qcdir/$name \
      --exclude $qcdir/$name.ac_gt_snps \
      --make-bed \
      --out $qcdir/$name.no_ac_gt_snps


#Prune Study Data --------------------------------------------------
plink --bfile  $qcdir/$name.no_ac_gt_snps \
      --indep-pairwise 50 5 0.2 \
      --out $qcdir/$name.no_ac_gt_snps

plink --bfile  $qcdir/$name.no_ac_gt_snps \
      --extract $qcdir/$name.no_ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/$name.pruned

#move log files
mv *.log $log


#Liftover study data back to hg37 -----------------------------------

# 1. Print UCSC BED file
awk '{print "chr" $1, $4 -1, $4, $2 }' $name.pruned.bim > ${name}_lift.bed

echo BED file DONE

# 2. Download (wget) chain files 38>37
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O hg38ToHg19.over.chain.gz
        
echo Chain file download DONE


# 3. liftBed: liftOver the bed file into build 37
liftOver ${name}_lift.bed hg38ToHg19.over.chain.gz ${name}_lift_out.bed ${name}_lift_out_unlifted.bed

echo Liftover DONE

# 4. Update data with liftover SNPs

# extract mapped variants snp IDs
awk '{print $4}' ${name}_lift_out.bed > ${name}_lift_map.snps
# extract updated positions
awk '{print $4, $3}' ${name}_lift_out.bed > ${name}_lift37coord.pos

echo Make SNP var and coord files DONE

#Update plink file with extracted snps and coordinates
plink2 --bfile ${name}.pruned --extract ${name}_lift_map.snps --update-map ${name}_lift37coord.pos --sort-vars --make-pgen --out ${name}_remap
#Error: Fixed-width .bed/.pgen output doesn't support sorting yet.  Generate a regular sorted .pgen first, and then reformat it
plink2 --pfile ${name}_remap --make-bed --out ${name}_remap

echo Update Plink files DONE

#Merge datasets so to keep snps that are just in 1000G data ---------------


#Generate Principle compenents ------------------------------------------

#Plot PCA in R ------------------------------------------------------------
