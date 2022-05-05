###Date: 2018-07-13
###Purpose: Plot user data projected on 1KG PCs, greying out 1KG individuals

setwd("~/RSV/data/post-imp/ancestry_38/qc_ancestry")

##Load packages
library(ggplot2)
library(colorspace)
library(readr)

##Initialise command line arguments
args <- commandArgs(TRUE)

##Load data
root <- args[1]

PCA <- read_delim("1KG_merged.eigenvec", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
PCA <- read_delim("1KG_merged.eigenvec", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
##Rename columns
colnames(PCA) <- c("FID","IID", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

#merge with Population data

Population <- read_delim("1KG_ID2Pop.txt")
merge <- full_join(PCA, pop, by = "IID")


#colour NAs as grey = our study data
#overlay that with the pop data
##Define colour palette
KG_Palette<-heat_hcl(length(unique(merge$Population)), h = c(300, 75), c. = c(35, 95), l = c(15, 90), power = c(0.8, 1.2), fixup = TRUE, gamma = NULL, alpha = 1)


##Print pairwise comparisons of PC1-5 to pdf
pdf("KG.pop_strat_PCA2.pdf")
with(merge, qplot(PC1,PC2,colour=Population) + scale_colour_manual(values = KG_Palette))
dev.off()

#merge %>% ggplot(aes(PC1,PC2, colour=Population)) + geom_point()

