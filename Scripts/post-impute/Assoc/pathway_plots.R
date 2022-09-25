# Pathway analysis plotting 20.09.22

# Set up ----------------------------------------------------------------------

setwd("~/RSV/data/post-imp/Assoc/nexus/no_sib")

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
library(forcats)
library(viridis)

path <- fread("pathway.csv")
path$parent_s <- as.factor(path$parent_s)

table(path$parent_s)


# Collapse similar partent factors ---------------------------------------------
# Recode signal transduction into one factor "Signal transduction"
# or coalesce table

signal <- levels(path$parent_s)[(str_detect(levels(path$parent_s), "Signal") == TRUE)]
path$parent_s <- 
  fct_collapse(path$parent_s, `Signal Transduction` = signal)

DNA <- levels(path$parent_s)[(str_detect(levels(path$parent_s), "DNA") == TRUE)]
path$parent_s <- 
  fct_collapse(path$parent_s, `DNA Replicarion` = DNA)

trans <- levels(path$parent_s)[(str_detect(levels(path$parent_s),
                                           regex("Transport", ignore_case = TRUE)) == TRUE)]

path$parent_s <- 
  fct_collapse(path$parent_s, `Vesicle mediated transport` = trans)


# Calculate weighted p-values -------------------------------------------------------

# mutate new var average p vale, grouping by parent factor
# The more obvs the more sig (so smaller pval) so this should roughly work
# Larger n will give a smaller mean p-val

path <-
  path %>%
  group_by(parent_s) %>%
  add_count(sort = TRUE) %>% 
  mutate(weighted_p = (mean(p_value) / n))


# Summarise -------------------------------------------------------------------------

summary <-
  path %>%
  group_by(parent_s) %>%
  summarise(mean_p = (mean(p_value)),
            count = n(),
            weighted_p = mean_p/count) %>%
  arrange(weighted_p)


  # I think p vals might add up indv p vals 


levels(summary$parent_s)

# Plot --------------------------------------------------------------------------------

# weighted p 
summary %>%
  ggplot(aes(x = fct_reorder(parent_s, weighted_p, .desc = TRUE), y = weighted_p, fill = weighted_p)) +
  geom_bar(stat = "identity") +
  xlab("\n") + ylab ("weighted p") + ggtitle("RSV Pathway Analysis") +
  coord_flip() +
  theme_minimal() + scale_fill_viridis(direction = -1)

# log weighted p

# add -log p val to df then plot
summary <- summary %>%
  mutate(log_p = -log2(weighted_p))
summary %>%
  ggplot(aes(x = fct_reorder(parent_s, log_p), y = log_p, fill = log_p)) +
  geom_bar(stat = "identity") +
  xlab("\n") +
  ggtitle("RSV Pathway Analysis") +
  ylab(expression(-log[2]~(italic(p)~weighted))) +
  coord_flip() +
  theme_minimal() + scale_fill_viridis(direction = 1) 

#i dont think this is a good way of plotting


# number of obvs
summary %>%
  ggplot(aes(x = fct_reorder(parent_s, mean_p), y = count, fill = mean_p)) +
  geom_bar(stat = "identity") +
  xlab("Parent Pathway\n") + ylab ("observations (n)") +
  coord_flip() +
  theme_minimal() + scale_fill_viridis()

summary %>%
  ggplot(aes(x = fct_reorder(parent_s, count), y = count, fill = mean_p)) +
  geom_bar(stat = "identity") +
  xlab("Parent Pathway\n") + ylab ("observations (n)") +
  coord_flip() +
  theme_minimal() + scale_fill_viridis(direction = -1)

# weighted p - count

summary %>%
  ggplot(aes(x = fct_reorder(parent_s, count), y = count, fill = weighted_p)) +
  geom_bar(stat = "identity") +
  xlab("Parent Pathway\n") + ylab ("observations (n)") +
  coord_flip() +
  theme_minimal() + scale_fill_viridis(direction = -1)



# Are number of observations counting unique snps-variant IDs or just subpathways?

# original pathway p is fishers exact test based on number of genes in top hits pathway, compared to the expected number of genes we would see in that pathway - so it is an enrichment p val, however doesnt take into account the strength of the GWAS p val from what I can tell


# Filter for sig pathways

sig <- filter(path, p_value <0.06)


sig_summary <-
  sig %>%
  group_by(parent_s) %>%
  summarise(mean_p = (mean(p_value)),
            count = n(),
            weighted_p = mean_p/count) %>%
  arrange(weighted_p)





# Just immune -------------------------------------------------------------------------

immune <- path[(path$parent_s %in% c("Immune System","Signal Transduction;Immune System","Disease")),]

immune <- geom_bar()


