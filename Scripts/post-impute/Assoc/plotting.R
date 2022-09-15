# Set up -----------------------------------------------------------------------

setwd("~/RSV/data/post-imp/Assoc")
plot_dir <- c("~/RSV/data/post-imp/Assoc/plots/")

top <-
  read.csv(file = "~/RSV/data/post-imp/Assoc/nexus/GWAS_summary_data.csv",
           stringsAsFactors = T)

  colnames(top)

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
  
  
# Calculate frequency of categories --------------------------------------------

levels(top$Predicted_Function) #Need to concat levels
levels(top$Feature_Type)
levels(top$Epigenome_Feature)
levels(top$Symbol) #131

freq <-
  setDT(top)[, .N, keyby = Feature_Type]

# Pie charts -------------------------------------------------------------------

bp <-
  ggplot(freq, aes(x = "", y = N, fill = Feature_Type)) +
  geom_bar(width = 1, stat = "identity")
pie <-
  bp + coord_polar("y", start = 0)

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )


his <-
  pie +
  blank_theme +
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#d645b2"))
his

pdf(file =(paste(plot_dir,"histone_pie.pdf",sep ="")))
his
dev.off()


# Manhattan Plot ---------------------------------------------------------------

gwas <- fread("RESV_total_score.glm.linear")

# rename cols for qqman
colnames(gwas) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1","TEST","OBS_CT","BETA",        
                    "SE","T_STAT","P","ERRCODE")

gwas$TEST <- as.factor(gwas$TEST)
gwas<- filter(gwas, TEST == "ADD") #Keep just genetic testing

write.csv(gwas, file = "RESV_gwas.csv")

#Plot Manhattan

pdf(file =(paste(plot_dir,"RESV_man.pdf",sep ="")))
manhattan(gwas, ylim = c(0, 8),
          cex = 0.6, cex.axis = 0.9,
          genomewideline = F,
          col = c("mediumseagreen", "darkolivegreen4"))
dev.off()



# Manhattan labels -------------------------------------------------------------

# Label top genes on manhatanna plot
# Can't use all labels as with LD there would be too many over laps
# So need to filter labels so there is only one, ideally the most sig gene


# Steps:
# 1) Arrange df so that it is grouped by factor var
# 2) slice_min(P) to take the lowest p value of each var
# 3) save as a new df, to be used as labels


test <- 
  top %>%
  group_by(Symbol) %>%
  slice_min(P)

test <- 
  top %>%
  group_by(Symbol) %>%
  summarise(minP = min(P))


#Public

public <-
  top %>%
  filter(!is.na(Public_GWAS_Trait))

del <- 
  top %>%
  filter(!is.na(Polyphen_Score))


# SNP freq plots ---------------------------------------------------------------

# Can't do bar charts as not case/control
# Plot SNP allele on x axis (group by)
# RESV score on y axis, maybe seperate violin plots



