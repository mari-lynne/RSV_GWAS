library(data.table)
library(dplyr)

setwd("/home/mari/RSV/data/post-imp/ancestry_38")

ref <- fread("all_hg38.psam")
ref_pheno <- fread("all_hg38.toMerge.fam")
join <- left_join(ref, ref_pheno, by ="IID", copy = TRUE)

#This is the fam we need to update "all_hg38.toMerge.fam"
("all_hg38.toMerge.fam")

#tables are in the same order so just keep the psam file rename to the fam file when writing table 
#just ad zeros column
ref <- ref %>% mutate(FID = rep(0,3202))
ref <- ref[,c(7,1:6)]
write.table(ref, file = "Ref_final.fam", row.names = F, quote=F, sep="\t")
#check then manually renamed to QC.fam

#Update sex ####
#Ref_final.fam also has sex pheno already
#Need to use pre-imp data as we filtered autosomes
library(data.table)
data <- fread("~/RSV/data/plink/rsv_QC2.fam") #this has IDs FID = IID hence how they got joined at somepoint

data <- fread("~/RSV/data/plink/rsv_QC2.fam")
covar <- fread("~/RSV/data/meta/covar.txt") #1 is male 2 is female
#fam file has .cel extension in ID name need to add string to update sex file 

sex <- select(covar,IID, Sex)
sex$IID <- paste(sex$IID,"CEL", sep = ".")

#add in fam colum duplicate IID column
sex$FID <- sex$IID
sex <- sex[,c(3,1,2)]


write.table(sex, file = "~/RSV/data/plink/update_sex.txt", quote = F, row.names = F, col.names = T, sep = "\t")


study <- fread("~/RSV/data/plink/study_samples.txt")
study$IID <- paste(study$IID, "CEL", sep = ".")

write.table(sex, file = "~/RSV/data/plink/study_samples.txt", quote = F, row.names = F, col.names = T, sep = "\t")


bim <- fread("~/RSV/data/plink/rsv_QC2.bim")
chr <- unique(bim$V1)
chr #chromosomes in bim file are labelled X, Y, XY, MT
table(bim$V1) #32495 X chromosomes, #1709 XY varients 8 Y varients

#--update-sex expects a file with FIDs and IIDs in the first two columns, and sex information (1 or M = male, 2 or F = female, 0 = missing

#-update-sex (Plink2) expects a file with sample IDs in front, and a sex information column.

#If there is a recognized header line (starting with '#FID' or '#IID'), it defaults to loading sex information from the first column titled 'SEX' (any capitalization); otherwise it assumes the 3rd column. To force a specific column number, use the 'col-num=' modifier.
#Only the first character in the sex column is processed. By default, '1'/'M'/'m' is interpreted as male, '2'/'F'/'f' is interpreted as female, and '0'/'N'/'U'/'u' is interpreted as unknown-sex. To change this to '0'/'M'/'m' = male, '1'/'F'/'f' = female, anything else other than '2' = unknown-sex, add the 'male0' modifier.


#Plots 2 make:
#Ancestry
#SNP function tally

#sex check plot

sexcheck <- fread("~/RSV/data/plink/rsv.sexcheck")

plot(sexcheck)
library(ggplot2)
ggplot(sexcheck, aes(x=F)) + geom_histogram(bins =75) + labs(title = "Plink Sex Check", x = " \nX-Chr Hetrozygozity (F-statistic)")




#?ylim  + ylim(-1.11,1.1
?geom_histogram
