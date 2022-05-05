library(data.table)
library(tidylog)
library(tidyr)
library(stringr)
library(dplyr)

#Post-impute plans:
#Files to make ####

- need to make covar file
- need to make pheno file
- filter indivuals in plink just to have case/control kids
- run pca/colour code by sex
- Run association analysis + Manhattan Plots 
- Redo 1000g PCA

#Check IDs in fam file 
setwd("~/RSV/data/post-imp/Assoc")

fam <- fread("RSV_imp_QC.fam")
#FIDs are 0's
#IIDs are A4634.cel file names (somewhere have been doubled, sep by _)

clin <- fread("~/RSV/data/meta/clinical_data_full_fin_comb.csv")
#full_IDs	Study.Subject.ID
#OX02-8860	OX020007 (clin data)


IDs <- fread("~/RSV/data/meta/oxford_sample_name_key.csv")
#link to biobank not array name

#Get recode file 
ID_link <- fread("~/RSV/data/meta/A4634_Sample_Annotation_and_QC_fails.csv")
#Other files are the other study data, need to speak to other sites about their metadata
#A-number (Same as fam), Sample_ID =OX01-8640-V1-DNA1
#This looks like either the biobank ID or Uniq_ID in other ID file
#Split visit into separate column

#Plan:
1) Update Clin file to inclue ID_link IDs via Sample ID > join with ID file (get doc ID)>
  2) Join by doc ID to clin file so it nnow has A-ID
    3) Fix fam file IDs (depduplicate), then update fam file directly in plink so we can match with clin data using A-ID

#Start with ID_link file   
A4634_0989
OX02-8497-DNA1
#Split into either biobank ID 8497 or docID OXO2-8497
#Need to basically trim last part
#I think code wise splitting Sample ID into all new columns via - delim is easiesy
#ID <- ID %>% mutate(Study = str_extract(META_ID, "\\w*(?=-)"))


ID_link <- ID_link %>% separate(Sample_ID, c("site", "Biobank", "Visit", "Sample"), sep = "-", remove=FALSE)

#Join with IDs file
IDs$Biobank <- as.character(IDs$Biobank)
test <- left_join(IDs, ID_link, by = "Biobank")
#Matched rows =153, #163 IDs leftover in ID file unjoined
library(tidylog)

#Lets try cat ID_link file site biobank to fullIDs to match clin data (Sep by -)
ID_link$fullIDs <- paste(ID_link$site, ID_link$Biobank, sep= "-")

clin2 <- left_join(clin, ID_link, by = "fullIDs")
#Matched rows =325 :))

#Add A-number ID to clin data file left join by Sample_ID
#Fix .fam file ####
fam2 <- str_split_fixed(fam$V2, "_", n = 2)
fam2 <- fam %>% mutate(Split_ID = str_extract(V2, "(?<=_).*"))
fam2 <- fam %>% mutate(Split_ID = sub("(.*?_.*?)_.*", "\\1", V2))
fam2 <- fam2 %>% mutate(Split_ID2 = str_replace_all(Split_ID, "\\.CEL", ""))

#reorder fam
fam2 <- fam2[,c(1,8,3:6)]
head(fam2)

#save.fam file
write.table(file = "~/RSV/data/post-imp/RSV_imp_QC.fam", fam2, row.names = F, col.names = T, quote=F, sep="\t")

#it is the psam file I want to modify
psam <- fread("RSV_imp_QC.psam")

psam2 <- select(fam2, Split_ID2, V5)
names(psam2) <- c("IID", "SEX")

write.table(file = "~/RSV/data/post-imp/RSV_imp_QC.psam", psam2, row.names = F, col.names = T, quote=F, sep="\t")

#Make pheno file ##########
#Reformat Clin2 columns
https://www.cog-genomics.org/plink/2.0/input#pheno
#Start with IID (plink2 doesnt need FID col)

clin2 <-clin2[, c(50:58,1:49)]
#clin2 <- dplyr::rename(clin2, IID = A-number)
names(clin2)[1] = "IID"
#just keep relevant cols
#clin3 <- clin2 %>% select(IID, Sex, Date.of.Birth, hospitalisation, respiratory_support, re)

#Do later for now just select --pheno-name RESV_total_score
clin3 <- filter(clin2, Comment != "FAIL QCCR")
#removed 81 rows
#clin3 <- drop_na(clin3, IID)

#Recode Sex ####
clin3$Sex <- clin3$Sex <- ifelse(clin3$Sex == "m", 1, ifelse(clin3$Sex == "f", 2, 0))
#Make exclusion list for QC fail (or actually just remove these from pheno file)
#remove na array numbers/FID
write.table(clin3, file = "pheno_all.txt", row.names = F, quote=F, sep="\t")
write.table(clin3, file = "covar_all.txt", row.names = F, quote=F, sep="\t")

#Keep list
https://www.cog-genomics.org/plink/2.0/filter
write.table(clin3$IID, file = "study_samples.txt", row.names = F, quote=F, sep="\t")

#Filtered COVAR/Pheno files simple
#Error: Invalid numeric token '03/06/2019' on line 6 of pheno_all.txt.

pheno <- clin3 %>% select(IID, RESV_total_score)
covar <- pca_covar %>% select(IID, PC1, PC2, PC3, PC4, Sex, has_pre_ex_con, baseline_age_at_visit)

write.table(pheno, file = "pheno.txt", row.names = F, quote=F, sep="\t")
write.table(covar, file = "covar.txt", row.names = F, quote=F, sep="\t")


#Other study sites ####
#Sites in clin - metadata
clin2$site <- as.factor(clin2$site)
levels(clin2$site)
summary(clin2$site)

#All sites in ID_link/which are all the samples genotyped
ID_link$site <- as.factor(ID_link$site)
levels(ID_link$site)
summary(ID_link$site)

#Double Check other folders for study site codes and metadata
#1210 samples in fam file across sites OX1,OX1A,OX2,SE01,SE01A,SE02,TU01A,UU01A,UU02,IM02,ED01A

#We have data from
- IM02
- OXO2
- SE02
- UU02

#Missing data from sites
- OXO1/A
- SE01/A
- TU01A
- UU01A
- ED01A

#Need to check if missing data can be found/if it is from the healthy control cohorts (but seems like a lot of healthy controls)
#Other things to look for - more clinical metadata such as pre-term, birth weight

clin2 <- rename(clin2, Geno_ID = IID)
write.csv(clin2, file = "current_metadata.csv", row.names = F)

#missing metadata
no_geno <- anti_join(clin, ID_link, by = "fullIDs")
#77 IDs in clin that were not genotyped (double check the casing)

no_meta <- anti_join(ID_link, clin, by = "fullIDs")
no_meta <- rename(no_meta, Geno_ID = `A-number`)
write.csv(no_meta, file = "missing_metadata.csv", row.names = F)


no_meta$site <- as.factor(no_meta$site)
levels(no_meta$site)
summary(no_meta$site)

#See email chain with Joe O2 samples = case/control, 01 samples are all from the infant active study (aka burden). Sample codes with 01A are from this arm of the study (children recruited at birth and sampled at all respiratory infections).
#So hopefully some of these kids get RSV but it is chance if they do or did



