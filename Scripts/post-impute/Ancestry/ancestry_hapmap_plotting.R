#Liftover and ancestry plotting:

#code based on https://meyer-lab-cshl.github.io/plinkQC/articles/HapMap.html and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025522/

#Aims: ####
#QC data for population structure
#check ancestry against hapmap data base
#can include pcas in association model as a covariate later

#Steps 
#Get HAPMAP data - use ftp function to download data from terminal
#update chromosome build of hapmap is in NCBI 36/hg18 we need it in 19 (GRCh37) using liftover
#merge with our datafile and compare ancestry
#along with HAPMAP data, need to QC our study data too, check for difficult to align snps, filter R2 before merging

#all steps apart from duplicates are done in bash/terminal 

#Liftover #####

#this tool updates the hapmap file so the build matches our GWAS data (I ended up using the online tool rather than code)
'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg18ToHg19.over.chain.gz'-O hg18ToHg19.over.chain.gz #don't need if using online tool https://genome.ucsc.edu/cgi-bin/hgLiftOver 

#download hapmap files, ftp = file transfer protocol
ftp=ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/
  prefix=hapmap3_r2_b36_fwd.consensus.qc.poly

#download files
wget $ftp/$prefix.map.bz2
bunzip2 $prefix.map.bz2

wget $ftp/$prefix.ped.bz2
bunzip2 $prefix.per.bz2

wget $ftp/relationships_w_pops_121708.txt

#I then renamed to hapmap3

#convert hapmap files to BED files for the liftover tool to work
#note this BED file is different to the plink .bed file (It took me way too long to realise this)
awk '{print "chr" $1, $4 -1, $4, $2 }' hapmap3.bim | \
sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > \
hapmap.tolift

#used online and got a .bed file with lifted co-ordinates NCBI build 36 to CGRCh37
# output reads Successfully converted 1439670 records: Conversion failed on 946 records. 
#download file and rename as HapMapIII_CGRCh37 #save in refdir folder

#get mapped positions
#awk is a command that you can use to deal with text or strings and prints them/reformats $ means the fourth field or column

# ectract mapped variants snp IDs
awk '{print $4}' HapMapIII_CGRCh37.bed > HapMapIII_CGRCh37.snps
# ectract updated positions
awk '{print $4, $3}' HapMapIII_CGRCh37.bed > HapMapIII_CGRCh37.pos


#update hapmap plink file snp coordinates using update map function
system("./plink --bfile hapmap3 \
--extract HapMapIII_CGRCh37.snps \
--update-map HapMapIII_CGRCh37.pos \
--make-bed \
--out hapmap3_remapped")

#output is --update-map: 1425733 values updated.
#Warning: Base-pair positions are now unsorted!
#Warning: 6287 het. haploid genotypes present (see hapmap3_remapped.hh ); many commands treat these as missing.

#copy and paste these files into GWAS directory for use
#cd to GWAS folder

#QC steps ####

#remove high LD snps in study data ####
#R2 > 0.2 to reduce computational complexity. Merging will struggle with regions in LD
#high-LD-regions.txt from #https://github.com/genepi-freiburg/gwas/tree/master/cleaning-pipeline/aux
#let me know if you have trouble downloading this!

system("./plink --bfile VAST_PATCH_QC3_rmdup --exclude high-LD-regions.txt --range --indep-pairwise 50 5 0.2 --out VAST_rmLD")
#creates VAST_rmLD.prune.in, the list of SNPs to be kept in the analysis. This step also removes SNPs from extended regions of high LD listed in high-LD-regions.txt.

#extract snps from my GWAS data
system("./plink --bfile VAST_PATCH_QC3_rmdup --extract VAST_rmLD.prune.in --make-bed --out VAST_PATCH.pruned")


#Create new bed file excluding SNPs in our VAST data which do not feature in the genotype data of the four original HapMap populations (HM3 data) saved as VAST.hapmap-snps
system("./plink --bfile VAST_PATCH.pruned --extract HapMapIII_CGRCh37.snps --make-bed --out VAST.hapmap-snps")

#remove duplicates in R #####
setwd('/home/mari/GWAS')
#install.packages("data.table")
library(data.table)
bim<-fread("VAST.hapmap-snps.bim")

bim.dups<-bim[duplicated(paste(bim$V1,bim$V4)),] 
duplicated_loci<-bim.dups$V2
write.table(file="duplicated_loci.txt",duplicated_loci,quote=FALSE,row.names = F,col.names = F)
system("./plink --bfile VAST.hapmap-snps --exclude duplicated_loci.txt --make-bed --out VAST_PATCH_rmdup ")#note had already previosuly removed duplicates, step not needed for me here

#Filter reference and study data for A->T or G->C SNPs ####
# awk script to find A→T and C→G SNPs. As these SNPs are more difficult to align and only a subset of SNPs is required for the analysis, we will remove them from both the reference and study data set. hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt downloaded from "https://github.com/genepi-freiburg/gwas/tree/master/single-pca
#rename to hapmap_no_atcg.txt

system("./plink --bfile VAST.hapmap-snps --extract hapmap_no_atcg.txt --make-bed --out VAST_cleaned")
#40129 pass QC

#do for hapmap files too
system("./plink --bfile hapmap3_remapped --extract hapmap_no_atcg.txt --make-bed --out hapmap_cleaned")
#1439670 variants loaded from .bim file. 973077 variants remaining.

#check chromosome mismatch ####
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
VAST_cleaned.bim hapmap_cleaned.bim | \
sed -n '/^[XY]/!p' > hapmap_cleaned.toUpdateChr

#theres only 3 snp mismatches 
./plink --bfile hapmap_cleaned \
--update-chr hapmap_cleaned.toUpdateChr 1 2 \
--make-bed \
--out hapmap_cleaned_chr

#mismatches continued 

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
VAST_cleaned.bim hapmap_cleaned.bim > hapmap_cleaned.toUpdatePos
#no mismatches detected :)

#Merging ####
# Merge the VAST files with the HapMap data
# This will give us a misssnip file that we can later use to check for mismatches/flips
system("./plink --bfile VAST_cleaned --bmerge hapmap_cleaned_chr --make-bed --out VAST_PATCH_merged")

#output: 40129 markers loaded from VAST_cleaned.bim.
#973077 markers to be merged from hapmap_cleaned.bim. Of these, 932948 are new, while 40129 are present in the base dataset.
#* If you believe this is due to strand inconsistency, try --flip with VAST_PATCH_merged-merge.missnp.

#flip mismatched snps in og datasets then remerge
system("./plink --bfile VAST_cleaned --flip VAST_PATCH_merged-merge.missnp --make-bed --out VAST_flipped.hapmap-snps")
#remerge with the flipped snps
system("./plink --bfile VAST_flipped.hapmap-snps --bmerge hapmap_cleaned_chr --make-bed --out VAST_PATCH_merged")

#only 1 variant left
#lets remove the bastard
system("./plink --bfile VAST_flipped.hapmap-snps --exclude-snp rs766173 --make-bed --out VAST_flipped.hapmap-snps")
#then remerge again as above
#Yay we have mergerd!!

#Plotting ####

#PCA on the merged data ####

./plink --bfile VAST_PATCH_merged --pca --out VAST_PATCH_merged

mv VAST_PATCH_merged $GWAS/plink_log

#PlinkQC ####
#R
setwd('/home/mari/GWAS')
library(plinkQC)
indir<-"/home/mari/GWAS"
qcdir<-"/home/mari/GWAS/QC"
name<-"VAST_cleaned"
path2plink <- "/home/mari/GWAS/plink"


refSamplesFile <- "/home/mari/GWAS/HapMap_ID2Pop.txt"
refColorsFile <- "/home/mari/GWAS/HapMap_PopColors.txt"
#downloaded from the author of the QC package's github https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/HapMap_ID2Pop.txt

exclude_ancestry <-
  evaluate_check_ancestry(indir=indir, name=name,
                          prefixMergedDataset="VAST_PATCH_merged",
                          refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
                                               sep=""), 
                          refColorsFile=paste(indir, "/HapMap_PopColors.txt",
                                              sep=""),
                          interactive=TRUE,
                          legend_text_size = 9,
                          legend_title_size = 11,
                          axis_text_size = 8,
                          axis_title_size = 9,
                          title_size = 12)

?check_ancestry


