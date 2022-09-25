# R set up ---------------------------------------------------------------------

setwd("~/RSV/data/transcriptomics")
plot_dir <-
  c("~/RSV/data/transcriptomics/plots") #where to save output plots

#load packages #
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(tidylog)
library(limma)
library(stringi)
library(janitor)
library(stringr)
library(splines)
library(edgeR)
library(BiocManager)
library(DESeq2)
library(Glimma)
library(RColorBrewer)
library(factoextra)
library(ggpubr)
library(MatrixEQTL)
library(biomaRt)



# Data set up ------------------------------------------------------------------

load(file = "expressionSetRma (03).Rda")

data = expressionSetRma
rm(expressionSetRma)

dims(data)
# Data = 573 samples

# Extract data from eset
samples <- sampleNames(data)
experiment <- experimentData(data)
# MIAME = minimum information about microarray experiment

exprs <- as.data.frame(data)
pData <- clean_names(data@phenoData@data)
fData <- clean_names(data@featureData@data)
pData <- setDT(pData, keep.rownames = TRUE)[]

# Filter expression data -------------------------------------------------------

# 1) Filter na genes
genes <- fData[(!is.na(fData$ensembl)),]
exprs <- exprs[,(colnames(exprs) %in% row.names(genes))]


# 1) Time point
table(pData$visit)

# V01 = 343
# V02 = 230
# 96 hours of symptom onset or 48 hours of admission to the hospital?
pData$visit <- str_replace_all(pData$visit, "\\d \\: ", "") # Tidy


v1_ids <- pData[(visit == "V01"),rn]
v1 <- pData[(visit == "V01"),]
v1_exprs <- exprs[(row.names(exprs) %in% v1_ids),]


# 2) Get snp data ---------- ---------------------------------------------------

assoc <- fread("~/RSV/data/post-imp/Assoc/RSV_Oxford_nosib.fam")

# DGE ID link file
annot <- read.csv(file = "Annotation3.csv")
annot <- clean_names(annot)
annot_v1 <- annot[(annot$name %in% v1_ids),]

# Geno ID file
# Need to get subject_ids for geno (currently in geno IDs)

geno_ids <- fread("~/RSV/data/meta/current_metadata.csv")
geno_ids <- dplyr::rename(geno_ids, subject_id = fullIDs)

geno_ids <- geno_ids[(geno_ids$Geno_ID %in% assoc$V2),] # Just get assoc samples from fam file


# Filter DGE data by assoc data -------------------------------------------------
# Filter dge annot file (link) by subject ID
# Then filter dge data

dge_geno <- annot_v1[(annot_v1$subject_id %in% geno_ids$subject_id),] # 277 samples in both

# Re - Filter exprs data and pheno
v1_exprs <- v1_exprs[(row.names(v1_exprs) %in% dge_geno$name),]
v1_pheno <- v1[(v1$rn %in% dge_geno$name),]

# Also have to remove/keep samples in plink that are in both data sets to eventually get genotyping info

geno_dge <- geno_ids[(geno_ids$subject_id %in% dge_geno$subject_id),] #293? must have duplicates
# remove duplicates
geno_dge <- geno_dge[(stri_duplicated(geno_dge$subject_id) == FALSE),]

# Keep these samples in plink
geno.keep <- assoc[(assoc$V2 %in% geno_dge$Geno_ID),2]
write.table(geno.keep, file = "~/RSV/data/meta/dge_v1_samples.txt", sep = "\t", quote = F, col.names = F, row.names = F)

# Filter for top SNPs
top <- fread("~/RSV/data/post-imp/Assoc/tophits22.txt")
snp.keep <- top$SNP
write.table(snp.keep, file = "~/RSV/data/transcriptomics/snp_keep.txt", sep = "\t", quote = F, col.names = F, row.names = F)

system("plink2 --bfile ~/RSV/data/post-imp/Assoc/RSV_Oxford_nosib --extract ~/RSV/data/transcriptomics/snp_keep.txt --keep ~/RSV/data/meta/dge_v1_samples.txt --make-bed --out ~/RSV/data/transcriptomics/RSV_dge")

save.image(file = "meqtl.RData")

# Prepare genotyping data/snp location data ------------------------------------

#Recode SNPs as 1's and 0's
# Open and run in terminal in ~/RSV/data/post-imp/Assoc/

# #Transpose table so SNPs are rows, participants are columns
# plink --bfile RSV_dge --allow-no-sex --recode oxford --out RSV_dge
# plink --data RSV_dge --allow-no-sex --recode A-transpose --out RSV_dge
# 
# # CHR	SNP	(C)M	POS	COUNTED	ALT
# cat RSV_dge.traw | cut -f2,7- > geno_file.txt

#cut relevant data from new file

# awk '{print $2, $1, $4}' RSV_dge.traw > snp_loc.txt


# Get gene location data -------------------------------------------------------

load(file = "meqtl.RData")
# Update array ids to gene names
colnames(v1_exprs) <- genes$symbol #from fData

# Download ensembl human gene data
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
# default build is 38

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes2 <- getBM(attributes=c('hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)

# Only assoc gene loc data with genes that are in our expression data
filtered_genes <- genes2[genes2$hgnc_symbol %in% (colnames(v1_exprs)),]

table(filtered_genes$chromosome_name)
filtered_genes <- filtered_genes[(filtered_genes$chromosome_name %in% c(1:22)),] # filter autosomes

# Write tables for meQTL -------------------------------------------------------
write.table(filtered_genes, file = "gene_loc.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Update Sample IDs to match geno IDs
colnames(v1_exprs) <- row.names(genes) # convert back as error if y columns aren't unique later

link <- left_join(dge_geno, geno_dge, by = "subject_id") %>%
  dplyr::rename(exprs_id = name) %>% dplyr::select(exprs_id, Geno_ID)

v1_exprs$exprs_id <- row.names(v1_exprs)

v1_exprs <- left_join(link, v1_exprs, by = "exprs_id")

row.names(v1_exprs) <-
  v1_exprs$Geno_ID #Update sample IDs to genotyping IDs
v1_exprs <- v1_exprs[, -c(1:2)] # Remove extra cols

# transpose exprs data
v1_exprs <- t(v1_exprs)

write.table(v1_exprs, file = "exprs.txt", sep = "\t", quote = F, col.names = T, row.names = T)

# Update covar data ------------------------------------------------------------

covar <- fread("~/RSV/data/post-imp/Assoc/pca_covar2.txt")
covar <- covar[(covar$IID %in% link$Geno_ID), ]

covar <- covar[, c(2:6,13:15)]

# Recode Sex
covar$Sex <- ifelse(covar$Sex == "Male", 1,
                    ifelse(covar$Sex == "Female", 2, NA))

covar <-  t(covar)
colnames(covar) <- covar[1,]
covar <- covar[-1,]
write.table(covar, file = "covar_eqtl.txt", sep = "\t", quote = F, col.names = T, row.names = T)


# Check order ------------------------------------------------------------------

# genotype = geno_file.txt
# expression = exprs.text
# covariates = covar_eqtl.txt
# gene location = gene_loc.txt
# SNP location = snp.loc.txt

# Check first 3 file column order

geno <- fread("~/RSV/data/transcriptomics/geno_file.txt")
covar <- fread("~/RSV/data/transcriptomics/covar_eqtl.txt")
exprs <- fread("exprs.txt")


# Clean geno names
colnames(geno) <- str_remove_all(colnames(geno), "0_")
geno <- dplyr::rename(geno, V1 = SNP)

# Order by geno
names.use <- names(geno)
exprs <- exprs[, ..names.use]
covar <- covar[, ..names.use]

write.table(geno, file = "eQTL/geno.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(exprs, file = "eQTL/exprs.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(covar, file = "eQTL/covar_eqtl.txt", sep = "\t", quote = F, col.names = T, row.names = F)

save.image(file = "geno.RData")


