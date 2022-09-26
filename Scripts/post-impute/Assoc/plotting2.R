# Plotting RSV 16/09/22 --------------------------------------------------------

setwd("~/RSV/data/post-imp/Assoc/nexus/no_sib")
load(file = "nexus.RData")
library(dplyr)
library(data.table)
library(ggplot2)
library(qqman)
library(ggrepel) 
library(tidyr)
library(tidylog)
library(stringi)
library(janitor)
library(stringr)
library(RColorBrewer)


# Feature pie charts -----------------------------------------------------------

main$feature_type[is.na(main$feature_type)] <- "N/A"
main$feature_type_class[is.na(main$feature_type_class)] <- "N/A"

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )


freq <-
  setDT(main)[, .N, keyby = feature_type_class]

# Pie charts -------------------------------------------------------------------

bp <-
  ggplot(freq, aes(x = "", y = N, fill = feature_type_class)) +
  geom_bar(width = 1, stat = "identity")
pie <-
  bp + coord_polar("y", start = 0)


his <-
  pie +
  blank_theme +
  scale_fill_brewer(palette = "Set3")
his

his1 <- pie +
  blank_theme +
  scale_fill_manual(values=c("#80B1D3", "lightgrey","#8DD3C7", "#FFFFB3"))
his1

(his1 + his) + plot_annotation(title = "Epigenome Features (SNPs)")

pdf(file = "epigenome.pdf", paper = "a4r")
(his1 + his) + plot_annotation(title = "Epigenome Features (SNPs)")
dev.off()


# Intron vs exon ----------------------

freq <-
  setDT(main)[, .N, keyby = annotation]

# New cats
# 3downstream,intronic
# 3downstream,non-coding intronic

main$annotation <- str_remove_all(main$annotation, "non-coding")
main$annotation <- str_replace_all(main$annotation, "codingintronicsyn", "intronic")
main$annotation <- str_replace_all(main$annotation, "intronicintronic", "intronic")
main$annotation <- str_replace_all(main$annotation, "3","3' ")
main$annotation <- str_replace_all(main$annotation, "5","5' ")

freq$annotation <- as.factor(freq$annotation)
levels(freq$annotation)
freq <- filter(freq, (levels(annotation) != ""))


bp <- ggplot(freq, aes(x = "", y = N, fill = annotation)) +
  geom_bar(width = 1, stat = "identity")

pie <-
  bp + coord_polar("y", start = 0)

snp_ft <-
  pie +
  blank_theme +
  scale_fill_manual(values = rev(brewer.pal(10, "Set3")))

snp_ft


#SNP genes
freq <-
  setDT(main)[, .N, keyby = type]

freq$type <- as.factor(freq$type)
levels(freq$type)
freq <- filter(freq, (levels(type) != "antisense"))
freq <- filter(freq, (levels(type) != "sense_overlapping"))

freq$type <- as.character(freq$type)

bp <- ggplot(freq, aes(x = "", y = N, fill = type)) +
  geom_bar(width = 1, stat = "identity")
bp

pie <-
  bp + coord_polar("y", start = 0)

snp_gene <-
  pie +
  blank_theme +
  scale_fill_manual(values = rev(brewer.pal(8, "Set3")))

snp_gene
snp_ft


rev(brewer.pal(8, "Set3"))

levels(freq$type) 
main$feature_type[is.na(main$feature_type)] <- "N/A"


levels(freq$type) <- gsub("^Non-coding$", "No-gene", levels(freq$type))

# Manhattan --------------------------------------------------------------------

# Steps:
# 1) Arrange df so that it is grouped by factor var
# 2) slice_min(P) to take the lowest p value of each var
# 3) save as a new df, to be used as labels

vignette("qqman")

DT <- data.table(main)

top <- DT[, .SD[which.min(p)], by = overlapped_gene]

top <- select(top,
              variation_id, chromosome, position, overlapped_gene, p)

# Manhattan labels -------------------------------------------------------------
assoc <- 
  fread("/home/mari/RSV/data/post-imp/Assoc/RESV_gwas.csv")


# make new gene col with genes and nearest genes (if no overlap)
# overwrite assoc$SNP with gene col
# Use annotateTOP = TRUE in man plot

save.image(file = "nexus.RData")


# Plot Manhattan ---------------------------------------------------------------

main$gene <- str_c(main$overlapped_gene, main$nearest_downstream_gene)
main$gene <- str_remove_all(main$gene, "None")

# Join to assoc table
genes <- select(main, chr, position, ref, alt, a1, gene)
genes <- rename(genes,
                CHR = chr,
                BP = position,
                REF = ref,
                ALT = alt,
                A1 = a1
                )

assoc2 <- left_join(assoc, genes, by = c("CHR","BP","REF","ALT","A1"))
#95 missing #1225 matching

assoc2$SNP <- assoc2$gene

assoc2 <- assoc2[,-1]

# Just need top value of each var

pdf(file = "RESV_man_genes.pdf", paper = "a4r")

manhattan(assoc2, snp = "SNP", ylim = c(0, 8),
          genomewideline = F,
          col = c("mediumseagreen", "darkolivegreen4"),
          annotatePval = 0.0001,
          main = "RSV Severity GWAS")

dev.off()


# Q-Q plot ---------------------------------------------------------------------

# Also check lamda/GIF 

qq2(assoc$P, main = "Q-Q plot of RSV GWAS p-values", col = "darkolivegreen4", xlim = c(0,7))

vignetes("qqman")



# should be in log file


# Prep Gtex data ---------------------------------------------------

gtex <- read.csv(file = "~/RSV/data/post-imp/Assoc/gtex.csv")

gtex <- drop_na(gtex) %>% filter(dbSNP != "none")

gtex <- gtex %>% filter(Symbol == "LSP1")

gtex$tissue <- rep("Lung", length(gtex$Symbol))

#gtex <- gtex[1:400,]

write.csv(gtex, file = "LSP1_lung.csv", row.names = F, quote = F)

top_main <- filter(main, p < 1e-5)

gtex <- gtex %>% filter(Symbol == "LSP1")
gtex$tissue <- rep("Whole_Blood", length(gtex$Symbol))
write.csv(gtex, file = "LSP1_blood.csv", row.names = F, quote = F)

gtex <- gtex %>% filter(Symbol == "MRC1")
gtex$tissue <- rep("Whole_Blood", length(gtex$Symbol))
write.csv(gtex, file = "MRC_blood.csv", row.names = F, quote = F)

# Plot Pathway analysis ---------------------------------------------------