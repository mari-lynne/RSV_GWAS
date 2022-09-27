# Tidy data ----------------------------------------------------------------

setwd("~/RSV/data/post-imp/Assoc/nexus/no_sib")

library(dplyr)
library(data.table)
library(ggplot2)
library(ggrepel) 
library(tidyr)
library(tidylog)
library(stringi)
library(janitor)
library(stringr)
library(RColorBrewer)
library(forcats)
library(viridis)


gwas_data <- clean_names(fread("gwas22.txt"))
gwas_data <- filter(gwas_data, !is.na(gwas_data$p))

# Reduce comp time by removing dense plots at the bottom which overlap

sig_data <- gwas_data %>% 
  subset(p < 0.01)
nsig_data <- gwas_data %>% 
  subset(p >= 0.01) %>%
  group_by(chr) %>% 
  sample_frac(0.05)
gwas_data <- bind_rows(sig_data, nsig_data)


# Calculate base pair coordinates -----------------------------
data_cum <- gwas_data %>% 
group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

gwas_data <- gwas_data %>%
select(-test,-obs_ct,-bp_add)

axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)
sig <- 1e-5

# Highlight genes --------------------------------------------

# Load data 
main <- read.csv(file = "GWAS_top_nosib.csv", header = T)
public <- clean_names(fread("gwas_RSV_nosib.txt"))
eqtl <- load(file = "~/RSV/data/transcriptomics/eQTL/eqtl_results.RData")
rm(snppos)


# Just remove dups for now
main <- main[(stri_duplicated(main$variation_id)==FALSE),]
main <- filter(main, !is.na(p))

# Update var_id format 
main$snp <- str_replace_all(main$variation_id, "\\/", "\\:") 
main$snp <- str_remove_all(main$snp, "(\\:1)$")
public$snp <- str_replace_all(public$variation_id, "\\/", "\\:") 
public$snp <- str_remove_all(public$snp, "(\\:1)$")


# variation_ids to highlight -------------------------------------------------------------
# filter so grouped by gene and then take the min val

DT <- data.table(main)
top <- DT[, .SD[which.min(p)], by = overlapped_gene]
top_snps <- top$snp

delete_snps <- main %>% subset(phred >= 12)
delete_snps <- delete_snps$snp
public_snps <- public$snp

cis_snps <- cis$snps
trans_snps <- trans$snps

# mutate table  -------------------------------------------------------------

# Add highlight and annotation information
gwas_data <- gwas_data %>%
  mutate(highlight_top = ifelse(snp %in% top_snps, "yes", "no")) %>%
  mutate(highlight_pub = ifelse(snp %in% public_snps, "yes", "no")) %>%
  mutate(highlight_del = ifelse(snp %in% delete_snps, "yes", "no")) %>%
  mutate(highlight_cis = ifelse(snp %in% cis_snps, "yes", "no")) %>%
  mutate(highlight_trans = ifelse(snp %in% trans_snps, "yes", "no"))

# Basic Plot ----------------------------------------------------

ggplot(gwas_data,
       aes(
         x = bp_cum,
         y = -log10(p),
         color = as_factor(chr),
         size = -log10(p)
       )) +
  geom_hline(yintercept = -log10(sig),
             color = "grey40",
             linetype = "dashed") +
  geom_point(alpha = 0.8, size = 1.3) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#44c17b", "#287349"), 22)) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = NULL,
       y = "-log(p)") + theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


# Plot with labelled eqtl genes ------------------------------------------------

# Get gene names
gene <- select(main, snp, overlapped_gene)
gwas_data <- left_join(gwas_data, gene, by = "snp")

eqtls <- rbind(cis, trans) %>% rename(snp = snps)
eqtls <- filter(eqtls, !stri_duplicated(gene))
eqtl_gene <- select(eqtls, snp, gene)
gwas_data <- left_join(gwas_data, eqtl_gene, by = "snp")

gwas_data <- filter(gwas_data, !is.na(gwas_data$p))


ggplot(gwas_data,
       aes(
         x = bp_cum,
         y = -log10(p),
         color = as_factor(chr),
         size = -log10(p)
       )) +
  geom_hline(yintercept = -log10(1e-4),
             color = "grey34",
             linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-5),
             color = "grey80",
             linetype = "dashed") +
  geom_point(alpha = 0.8, size = 1.3) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#44c17b", "#287349"), 22)) +
  scale_size_continuous(range = c(0.5, 3)) + geom_label_repel(
    data = subset(gwas_data, highlight_trans == "yes"),
    aes(label = gene),
    size = 2.5, color = "#AD3A8A",
  max.overlaps = 20, nudge_y = 0.005) +
  geom_label_repel(
    data = subset(gwas_data, highlight_cis == "yes"),
    aes(label = gene),
    size = 2.5, color = "#D28423",nudge_y = 0.1, nudge_x = 0.05
  ) +
  labs(x = NULL,
       y = "-log(p)") + theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# Anotate ----------------

main_all <- read.csv(file = "GWAS_top_nosib.csv", header = T)
main_all$snp <- str_replace_all(main_all$variation_id, "\\/", "\\:") 
main_all$snp <- str_remove_all(main_all$snp, "(\\:1)$")

eqtl_map <- main_all[main_all$snp %in% eqtls$snp, ]

rm(main_all)


save.image(file = "manhattan.RData")

# Public hits ---------------------------------------------
rm(gwas_data)

gene <- select(public, snp, genes,trait) %>% filter(genes != "None")

gene <- mutate(gene, gene_trait = paste(gene$genes, gene$trait, sep = " - "))

# delete extra lisp for now
gene <- gene[-c(8:9,12),]

gwas_data <- left_join(gwas_data, gene, by = "snp")
gwas_data <- filter(gwas_data, !is.na(gwas_data$p))


# Plot
ggplot(gwas_data,
       aes(
         x = bp_cum,
         y = -log10(p),
         color = as_factor(chr),
         size = -log10(p)
       )) +
  geom_hline(yintercept = -log10(1e-4),
             color = "grey34",
             linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-5),
             color = "grey80",
             linetype = "dashed") +
  geom_point(alpha = 0.8, size = 1.3) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#44c17b", "#287349"), 22)) +
  scale_size_continuous(range = c(0.5, 3)) +
  geom_label_repel(
    data = subset(gwas_data, highlight_pub == "yes"),
    aes(label = gene_trait),
    size = 2.5, color = "sienna3",nudge_y = 0.6, nudge_x = 0.5, direction = "both"
  ) +
  labs(x = NULL,
       y = "-log10(p)", title = "Public GWAS Hits") + theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )






