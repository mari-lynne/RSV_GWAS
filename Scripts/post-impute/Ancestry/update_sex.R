library(data.table)
library(dplyr)
library(ggplot2)

setwd("/home/mari/RSV/data/plink")

#Update sex ####
#Ref_final.fam also has sex pheno already
#Need to use pre-imp data as we filtered autosomes

#Check covar file ####
data <- fread("rsv_QC2.fam") #this has IDs FID = IID hence how they got joined at somepoint
covar <- fread("~/RSV/data/meta/covar.txt") #1 is male 2 is female

#Make update sex file ####
#fam file has .cel extension in ID name
sex <- select(covar,IID, Sex)
#Add .CEL string to update sex file using paste
sex$IID <- paste(sex$IID,"CEL", sep = ".")

#add in fam column by duplicating IID column
sex$FID <- sex$IID
#Reorer so is FID, IID, Sex (as requested by plink)
sex <- sex[,c(3,1,2)]
write.table(sex, file = "~/RSV/data/plink/update_sex.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#--update-sex expects a file with FIDs and IIDs in the first two columns, and sex information (1 or M = male, 2 or F = female, 0 = missing
#-update-sex (Plink2) expects a file with sample IDs in front, and a sex information column.
#If there is a recognized header line (starting with '#FID' or '#IID'), it defaults to loading sex information from the first column titled 'SEX' (any capitalization); otherwise it assumes the 3rd column. To force a specific column number, use the 'col-num=' modifier.
#Only the first character in the sex column is processed. By default, '1'/'M'/'m' is interpreted as male, '2'/'F'/'f' is interpreted as female, and '0'/'N'/'U'/'u' is interpreted as unknown-sex. To change this to '0'/'M'/'m' = male, '1'/'F'/'f' = female, anything else other than '2' = unknown-sex, add the 'male0' modifier.


#Remake study sample keep list ####
#Remake study samples keep list to include .cel in ID also
study <- fread("~/RSV/data/plink/study_samples.txt")
study$IID <- paste(study$IID, "CEL", sep = ".")
write.table(sex, file = "~/RSV/data/plink/study_samples.txt", quote = F, row.names = F, col.names = T, sep = "\t")


#Check chromosome annotation in bim ####
#rsv_QC2 results:
#Chromosomes labelled X, Y, XY, MT
#32495 X chromosomes, #1709 XY variant 8 Y variant

#Redo with LD pruning: rsv_QC2.bim > LD > rsv_sex_clean.bim
bim <- fread("~/RSV/data/plink/rsv_sex_clean.bim") 
#Chromosomes in bim file are labelled 23, 24, 25, 26
table(bim$V1) #7677 X chromosomes, #979 XY variant, 1 Y variants, 25 MT variants

#sex check plot ############################

sexcheck <- fread("~/RSV/data/plink/rsv.sexcheck")
ggplot(sexcheck, aes(x=F)) + geom_histogram(bins =25) + labs(title = "Plink Sex Check -hh", x = " \nX-Chr Homozygosity (F-statistic)") #600W, H=385
#below 0.2 = Female, above 0.8 = male (looks right from the data)

#Heterozygosity check #####################
#145 sample mis-sexed seems like a lot!
#136021 het. haploid genotypes present (see rsv_sex_clean.hh ); many commands treat these as missing.

hh <- fread("rsv_sex_clean.hh", header = F)
#Family ID, Within-family ID,Variant ID

bim <- rename(bim, SNP = V2)
hh <- rename(hh, SNP = V3)

chr_check <- left_join(hh, bim, by = "SNP")
table(chr_check$V1.y)
#Chr 23 = X chromosome where all the hetezygous snps are, isn't this what we would expect?

#Remove with --set-hh-missing
#Now 65 mismatches remain
#Redo plot > it seems it correctly assigned more females (F >0.8)
#Remaining problems are all ped-sex = 2 (F), snp sex = 1 (M) 
#So females XX being incorrectly assigned X-Y, they have low homozygosity

#If heterozygous haploid calls still remain, the most likely cause is nonmissing female genotype calls on the Y chromosome; others have reported that this is fairly common.  A quick way to check the number of these is to just load the Y chromosome with e.g. "plink --bfile semi_clean_fileset --chr 24 --freq".  If all the heterozygous haploid errors are on the Y chromosome, you can safely clobber them with --make-bed + --set-hh-missing.  (If some are on the X, --set-hh-missing *might* still be okay, but I'd need to know more about the data source and the --check-sex report to be sure.)

#Check pheno file ####

fail <- filter(sexcheck, STATUS == "PROBLEM")


meta <- fread("~/RSV/data/meta/covar_all.txt")
meta_og <- fread("~/RSV/data/meta/current_metadata.csv")

meta$IID <- paste(meta$IID, "CEL", sep = ".")
meta_og$IID <- paste(meta_og$Geno_ID, "CEL", sep = ".")

check <- left_join(fail, meta_og, by = "IID")

#Double check ID matching worked
clin <- fread("~/RSV/data/meta/clinical_data.csv")
IDs <- fread("~/RSV/data/meta/oxford_sample_name_key.csv")
ID_link <- fread("~/RSV/data/meta/A4634_Sample_Annotation_and_QC_fails.csv")

#manual check problem samples
#From the top of fail

grep("A4634_0169", ID_link$A-number)

#A4634_0169 SE02-0013-V01-RB1 0.98 99.6 #sex is f in clin file
#x4 manual checks and all sex are the same in original metadata = f so not a recoding problem




#Double check variant quality (should have been qc'd)
#Redo steps
#Remove XY region - Just select XY chrs
#Add QC thresholds
#LD prune
#Then run split-x 37
#Then remove (set missing hh)
#Then run sex-check with y counts

hh <- fread("rsv_redo.hh", header = F)
#Family ID, Within-family ID,Variant ID
bim <- rename(bim, SNP = V2)
hh <- rename(hh, SNP = V3)
chr_check <- left_join(hh, bim, by = "SNP")
table(chr_check$V1.y)
#Chr 23 = X chromosome where all the hetezygous snps are, isn't this what we would expect?


#More file cleaning
#Load SNP info
bim <-fread("rsv_QC2.bim")
#Filter SNPs in autosomal chromosomes")
xy <- bim %>% filter(V1 == "Y" | V1 == "XY" | V1 == "X")

#Remove duplicates 
duplicated_snps <- xy[duplicated(paste(xy$V1,xy$V4)),]
#remove snps where A > T, C > G
#PLINK can't tell what strand to merge/read these SNPs from
atcg_snps <-  xy[((xy$V5 == "A") & (xy$V6 == "T")) |
                    ((xy$V5 == "T") & (xy$V6 == "A")) |
                    ((xy$V5 == "C") & (xy$V6 == "G")) |
                    ((xy$V5 == "G") & (xy$V6 == "C"))]
#Exclusion list
exclude <- rbind(duplicated_snps, atcg_snps)
exclude <- exclude$V2
write.table(file="exclude_snps_xy.txt", exclude, sep = "\t", quote=F, row.names = F, col.names = F)

#No joy - biostars may be an affymetrix array issue
#https://www.biostars.org/p/207206/#9521902


