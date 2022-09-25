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

pData <- clean_names(data@phenoData@data)
fData <- clean_names(data@featureData@data)
exprs <- as.data.frame(data@assayData$exprs)


# GWAS data set up  --------------------------------------------------------------------

# Assoc log file:
# 314 samples (RESV cases) RSV_Oxford_nosib
# Covar/metadata in pca_covar2.txt

gwas <- # read.csv(file = "~/RSV/data/meta/current_metadata.csv")
gwas <- gwas %>% drop_na(Geno_ID) # Filtering samples that don't have a genotyping ID
gwas <- clean_names(gwas)



# START DGE data set-up --------------------------------------------------------------

exprs <- as.data.frame(data)
pData <- clean_names(data@phenoData@data)
fData <- clean_names(data@featureData@data)
pData <- setDT(pData, keep.rownames = TRUE)[]

# pData has visit data, but no other meta data, get by left joining to GWAS
# Need to add IDs

# pData <- cbind(pData, samples$row.names.samples) ?

gwas <-
  dplyr::rename(gwas, subject_id = full_i_ds)

gwas <-
  left_join(gwas, annot, by = "subject_id")

colnames(gwas) <-
  str_replace_all(colnames(gwas), ".x", "")

gwas <- gwas[,-c(6:9)]
gwas <- gwas[,-c(57:61)]
colnames(pData)

gwas <- 
  dplyr::rename(gwas, exp_id = name)
pData <-
  dplyr::rename(pData, exp_id = V2)
pData <- 
  pData[,c(10,1:9)]

test <-
  left_join(pData, gwas, by = "exp_id")

strict_left_join <- function(x, y, by = NULL, ...) {
  by <- common_by(by, x, y)
  if (any(duplicated(y[by$y]))) {
    stop("Duplicate values in foreign key")
  } else
    left_join(x, y, by = by, ...)
}


# Check for duplicate names ------------------------------------------

g2 <- gwas %>% mutate(dup_check = stri_duplicated(gwas$subject_id))


g2 <- gwas %>% mutate(dup_check = stri_duplicated(gwas$exp_id))
#highlights first duplicate entry
#remove duplicates
test <- g2 %>% filter(dup_check == TRUE)

g2 <- g2 %>% filter(dup_check == TRUE)
test <-
  left_join(pData, g2, by = "exp_id")

# finish later 

#Data vis -------------------------------------------------------

#Time point MDS
group = as.factor(data$samples$time)
col.group <- group
#Set colours
palette_Dark2 <- colorRampPalette(brewer.pal(8, "Set2"))
levels(col.group) <- palette_Dark2(length(unique(data$samples$time)))
col.group <- as.character(col.group)
#Plot MDS
pdf(paste(plot_dir, "time_mds.pdf", sep = ""))
#starts writing a PDF to file
plotMDS(lcpm, labels = group, col = col.group)
title(main = "Time Point") #Some clustering around TD already :)
dev.off()

