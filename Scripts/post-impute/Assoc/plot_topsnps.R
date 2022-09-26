# Plot top snps:

# Aims:
# Visualise top snps by genotype frequency
# SPlit by case control + across RESV score 

# Steps:
# Use recoded oxford file for 01 data
# Merge with pheno data by geno ID
# Plot

# Set up -----------------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(tidylog)
library(stringi)
library(janitor)
library(stringr)
library(splines)
library(RColorBrewer)
library(ggpubr)
library(biomaRt)


setwd("~/RSV/data/transcriptomics")

geno <- fread("~/RSV/data/transcriptomics/geno_file.txt")
pheno <- fread("~/RSV/data/post-imp/Assoc/pheno2.txt")
covar <- fread("~/RSV/data/transcriptomics/eqtl_covar.txt")



# Format data ------------------------------------------------------------------
geno <- as.data.frame(t(geno))
row.names(geno) <- str_remove_all(row.names(geno), "0_")
colnames(geno) <- geno[1,]
geno <- geno[-1,]

geno$IID <- row.names(geno)
geno <- left_join(geno, pheno, by = "IID")

# Plot -------------------------------------------------------------------------

# SNPs of interest 
#chr15:54797788:T:G

geno %>%
  ggplot(aes(y = RESV_total_score, x = `chr15:54797788:T:G`)) +
  geom_violin() + theme_pubclean() + scale_colour_brewer(palette="Dark2")



ggviolin(geno, x = "chr15:54797788:T:G",
         y = "RESV_total_score",
         color = "chr15:54797788:T:G", palette = "jco",
         ylab = "RESV Score", 
         add = "median_iqr")


ggviolin(geno, x = "chr15:54797788:T:G",
         y = "RESV_total_score",
         color = "chr15:54798159:C:T", palette = "jco",
         ylab = "RESV Score", 
         add = "median_iqr")



# Loop plots
# list of eqtl snps
