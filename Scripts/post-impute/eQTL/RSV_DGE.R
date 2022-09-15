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
# Data1 = 1107 samples

# Extract data from eset
samples <- sampleNames(data)
experiment <- experimentData(data)
# MIAME = minimum information about microarray experiment

pData <- clean_names(data@phenoData@data)
fData <- clean_names(data@featureData@data)
exprs <- as.data.frame(data@assayData$exprs)


# Check IDs --------------------------------------------------------------------
gwas <- read.csv(file = "~/RSV/data/meta/current_metadata.csv")

pData <- setDT(pData, keep.rownames = TRUE)[]
gwas <- gwas %>% drop_na(Geno_ID)
gwas <- clean_names(gwas)


# Check annotation file --------------------------------------------------------

annot <- read.csv(file = "Annotation3.csv")
annot <- clean_names(annot)


gwas <-
  dplyr::rename(gwas, subject_id = full_i_ds)


test <- left_join(gwas, annot, by = "subject_id")

#Do need to fix, resceu ID, subject_ID has array names mid way thru
write.table(annot, file = "samples.txt", row.names = F)
system(awk '{print $2}' samples.txt |cut -d "-" -f1,2 |sed -e 's/"//g'> fullIDs.txt)

IDs <- read.csv(file = "fullIDs.txt")

names(IDs) <- c("subject_id")

annot <- annot[,-2]
annot <- cbind(annot,IDs)
rm(IDs)

annot <- left_join(gwas,annot, by = "subject_id")
#rows only in x 97 > rows only in y(  1), > matched rows     521 
rmdp(test)
#filter data after DGE analysis

# DGE data set-up --------------------------------------------------------------

exprs <- as.data.frame(data)
pData <- clean_names(data@phenoData@data)
fData <- clean_names(data@featureData@data)

# pData has visit data, but no other meta data, get by left joining to GWAS
# Need to add IDs

pData <- 
  cbind(pData, samples$row.names.samples.)

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

