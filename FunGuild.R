library(readr)
library(dplyr)
library(tidyr)
# Thilina's project
df <- read_csv("~/Downloads/FUNGuild_results_Group.csv")
#otu <- df[,3:70]
otu <- df[,3:13]
#sums <- colSums(otu)
# install.packages('MultIS')
library('MultIS')
otu1 <- convert_columnwise_relative(as.matrix(otu))
df[,3:70] <- otu1
df$otu_sum <- rowSums(df[,3:70])
df <- df %>% select(ASV_ID, growthForm, otu_sum)
df.1 <- tidyr::spread(df, growthForm, otu_sum, fill = 0)
write.csv(df.1, 'FUNGuild_relative_abundance.csv')
# bar plot for growthform
gf <- as.data.frame(otu1)
gf$GF <- df$growthForm
gf_long <- reshape2::melt(gf, id.var = "GF")
gf_long.1 <- gf_long %>% group_by(GF, variable) %>% summarise(Sum_value = sum(value))
library(ggplot2)
ggplot(gf_long.1, aes(fill=GF, y=Sum_value, x=variable)) + 
  geom_bar(position="stack", stat="identity")
gf_long.1 %>% filter(GF == 'Coralloid')
# make into a wide form for ordination, using spread.

# correlation analysis
# https://microbiome.github.io/tutorials/Heatmap.html
# In the example, the two matrix are :
# x = asv * sample
# y = lipids * sample
# The result was: asv * lipids
# In our study, Thilina wants
# x sample * asv    

# Ficus project
guilds <- read_delim("Downloads/guilds.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)


# command: grep -i arbuscular guilds.txt > amf.guilds.txt
AMF_guilds <- read_delim("Downloads/AMF.guilds.txt", 
                              delim = "\t", escape_double = FALSE, 
                              col_names = FALSE, trim_ws = TRUE)
colnames(AMF_guilds) <- colnames(guilds)
# column: guilds
# ASV #: 73
# intersection: no pathogens
AMF_guilds <- guilds %>% filter(guild == 'Arbuscular Mycorrhizal')
# 2. Plant pathogens
# key words: "plant pathogen"
# command: grep -i 'plant pathogen' guilds.txt > plant_pathogen.guilds.txt
pp_guilds <- read_delim("Downloads/plant_pathogen.guilds.txt", 
                                         delim = "\t", escape_double = FALSE, 
                                         col_names = FALSE, trim_ws = TRUE)
colnames(pp_guilds) <- colnames(guilds)
pp_guilds$trophicMode <- as.factor(pp_guilds$trophicMode)
summary(pp_guilds$trophicMode)
# overlaps a lot with endophytes
guilds$trait <- as.factor(guilds$trait)
summary(guilds$trait)
pp_all <- guilds %>% filter(grepl('Plant Pathogen', guild, ignore.case=FALSE))
# 3. endophytes
endo_guilds <- guilds %>% filter(grepl('Endophyte', guild, ignore.case=FALSE))

