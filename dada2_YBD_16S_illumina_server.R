library(dada2)
library(ggplot2)
library(phyloseq)
library(phangorn)
library(DECIPHER)
library(dplyr)

library(readr)
library(vegan)
library(ggrepel)
library(tibble)
library(devtools)
source('https://raw.githubusercontent.com/vqv/ggbiplot/master/R/ggbiplot.r')

source('https://raw.githubusercontent.com/vqv/ggbiplot/master/R/ggscreeplot.r')

library(FactoMineR)
library(factoextra)
library(ggsci)
path <- '/media/hfan/T5/YBD/2019Sep/metabarcode/16S/fastq/'
files = list.files(path, pattern = '.fq', full.names = TRUE)
# plotQualityProfile(files)
# we also need the sample names
files_list = list.files(path, pattern = '.fq')
sample_names <- sapply(strsplit(files_list,"\\."),
                       `[`,  # backtick (or back quote, `), [ for exact match while [[ can extend
                       1) # indexing by character, extracts the first element of a subset
sample_data <- read.csv('/media/hfan/T5/YBD/2019Sep/YBD_metadata.csv',
                        header=TRUE)

rownames(sample_data) <- sample_data$Sample
filtered_path <- '/media/hfan/T5/YBD/2019Sep/metabarcode/16S/trimmed'
trim_files <- file.path(filtered_path,
                              paste0('trim_', sample_names, ".fastq.gz"))



# Primer used: 515F-806R, 292bp, V4 of 16S
# need to filter the long&crappy reads (using python), later
out <- filterAndTrim(files, trim_files, maxLen = 300, maxN=0,
                     maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE,
                     multithread=TRUE, trimLeft = 23, trimRight = 20)
# 
errors <- learnErrors(trim_files, multithread=TRUE)
saveRDS(errors,'/media/hfan/T5/YBD/2019Sep/metabarcode/16S/dada2_errors.rds')
# time demanding. But why only 12 samples

derep <- derepFastq(trim_files, verbose=TRUE)

# name the derep-class objects by the sample names
names(derep) <- sample_names

dada_1 <- dada(derep, err=errors, multithread=TRUE)
seq_table <- makeSequenceTable(dada_1)
dim(seq_table)

# inspect distribution of sequence lengths
table(nchar(getSequences(seq_table)))
seqtab2 <- seq_table[,nchar(colnames(seq_table)) %in% 252:253]
seq_table_nochim <- removeBimeraDenovo(seqtab2, method='consensus',
                                       multithread=TRUE, verbose=TRUE)
dim(seq_table_nochim)

# which percentage of our reads did we keep?
sum(seq_table_nochim) / sum(seq_table)

get_n <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(trim_files, get_n), sapply(dada_1, get_n),
               rowSums(seq_table), rowSums(seq_table_nochim))

colnames(track) <- c('input', 'filtered', 'denoised', 'merged', 'tabled',
                     'nonchim')
rownames(track) <- sample_names
head(track)

taxa <- assignTaxonomy(seq_table_nochim,
                       '/media/hfan/T5/YBD/silva_nr_v128_train_set.fa.gz',
                       multithread=TRUE)

# time consuming

#seq_table <- readRDS('seq_table.rds')
taxa <- addSpecies(taxa, '/media/hfan/T5/YBD/silva_species_assignment_v128.fa.gz')
saveRDS(taxa, '/media/hfan/T5/YBD/2019Sep/metabarcode/16S/dada2_taxa.rds')
#saveRDS(seq_table, 'seq_table.rds')
#saveRDS(seq_table_nochim, 'seq_table_nochim.rds')
#taxa <- readRDS('taxa.rds')
taxa_print <- taxa  # removing sequence rownames for display only
rownames(taxa_print) <- NULL
head(taxa_print)

sequences <- getSequences(seqtab2)
names(sequences) <- sequences  # this propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA) # time consuming
phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
saveRDS(phang_align,'/media/hfan/T5/YBD/2019Sep/metabarcode/16S/dada2_phang_align.rds')
#dm <- dist.ml(phang_align)
#treeNJ <- NJ(dm)  # note, tip order != sequence order
# takes a long time! Also I don't need this thing.
#saveRDS(treeNJ, 'treeNJ_half.rds')
#treeNJ <- readRDS('treeNJ.rds')
#phang_align <- readRDS('phang_align.rds')
# pml computes the likelihood of a phylogenetic tree given a sequence alignment and a model. optim.pml optimizes the different model parameters.
#fit = pml(treeNJ, data=phang_align)

## negative edges length changed to 0!

#fitGTR <- update(fit, k=4, inv=0.2)
#saveRDS(fitGTR, 'fitGTR_half.rds')
#fitGTR <- readRDS('fitGTR.rds')
# construct ML tree
# fitGTR <- optim.pml(fitGTR, model='GTR', optInv=TRUE, optGamma=TRUE,
#                    rearrangement = 'stochastic',
#                    control = pml.control(trace = 0))
#Error: vector memory exhausted (limit reached?)
# First, check whether you are using 64-bit R. 
# I gave up! Too many tips.
detach('package:phangorn', unload=TRUE)

#seq_table_nochim <- readRDS('seq_table_nochim.rds')
# The blow error message was caused by IDing with integers
### Error in validObject(.Object) : invalid class “phyloseq” object: 
### Component sample names do not match.
### Try sample_names()

physeq <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE),
                   sample_data(sample_data),
                   tax_table(taxa))

physeq1 <- physeq
OTU1 <- as(otu_table(physeq1), "matrix")
OTU2 <- t(OTU1)
taxa_names(physeq1) <- paste0("SV", seq(ntaxa(physeq1)))
OTU3 <- as(otu_table(physeq1), "matrix")
OTU4 <- as.data.frame(t(OTU3))
OTU4$OTU <- rownames(OTU4)
OTU4 <- OTU4[1:(length(OTU4)-1)]
write.table(OTU4,'/media/hfan/T5/YBD/2019Sep/metabarcode/16S/OTU4.csv',row.names = FALSE, sep = ',', quote = FALSE)
# save taxa to file
seq_taxa <- merge(taxa, OTU2, by=0)
seq_taxa <- seq_taxa %>% rename(ASV = Row.names)
write.table(seq_taxa, '/media/hfan/T5/YBD/2019Sep/metabarcode/16S/seq_taxa.csv',
            row.names = FALSE, sep = ',', quote = FALSE)
# save relative abundance
OTU1_relative <- as(otu_table(physeq_relative), "matrix")
OTU2_relative <- t(OTU1_relative)
seq_taxa <- merge(taxa, OTU2_relative, by=0)
seq_taxa <- seq_taxa %>% rename(ASV = Row.names)
write.table(seq_taxa, '/media/hfan/T5/YBD/2019Sep/metabarcode/16S/seq_taxa_relative_abundance.csv',
            row.names = FALSE, sep = ',', quote = FALSE)

taxa_names(physeq1) <- paste0(">SV", seq(ntaxa(physeq1)))
OTU5 <- as(otu_table(physeq1), "matrix")
OTU6 <- as.data.frame(t(OTU5))
OTU6$Seq <- colnames(OTU1)
OTU6 <- OTU6 %>% select(Seq)
write.table(OTU6, file = '/media/hfan/T5/YBD/2019Sep/metabarcode/16S/OTU6.fa', quote = FALSE, sep = '\n', col.names = FALSE)


plot_richness(physeq, x='Species', measures=c('Shannon', 'Fisher'), color='Habitat') +
    theme_minimal()

ord <- ordinate(physeq, 'MDS', 'euclidean')
plot_ordination(physeq, ord, color='Species', shape = 'Habitat', 
                title='PCA of 16S data') +
    theme_minimal()

ord <- ordinate(physeq, 'NMDS', 'bray')
plot_ordination(physeq, ord, type='samples', color='Species', shape = 'Habitat',
                title='PCA of the samples from the MiSeq SOP') +
    theme_minimal()

top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]
physeq_relative <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
physeq_top20 <- prune_taxa(top20, physeq_relative)

plot_bar(physeq_top20, x="Sample", fill='Genus') +
    facet_wrap(~Habitat, scales='free_x') +
    theme_minimal()

plot_bar(physeq_top20, x="Sample", fill='Family') +
    facet_wrap(~Species, scales='free_x') +
    theme_minimal()

# HCPC for microbial community
# microbial
# lots of NAs at different level, even Order! any(is.na(seq_taxa$Order)) ==  TRUE
# So what do we do? Rename them to "Unclassified".
level <- 'Family'
data_F <- seq_taxa 
data_F$Family[is.na(data_F$Family)] <- "Unclassified" #need to aggregate unclassified ones
data_F <- data_F %>% select(Family, Cow_Riv_2:Goat_Thi_3)
# aggregate by family
data_F_sum <- as.data.frame(data_F %>% group_by(Family) %>% summarise_all(sum))
data_G <- seq_taxa
data_G$Genus[is.na(data_G$Genus)] <- "Unclassified"
data_G <- data_G %>% select(Genus, Cow_Riv_2:Goat_Thi_3)
data_G_sum <- data_G %>% group_by(Genus) %>% summarise_all(sum)

# NMDS and HCPC
rownames(data_F_sum) <- data_F_sum$Family
data_F_sum <- data_F_sum %>% select(-Family)
# HCPC
N_matrix <- t(as.matrix(data_F_sum))
res.pca <- PCA(N_matrix,ncp=3, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)
fviz_dend(res.hcpc, 
          cex = 0.9,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 35      # Augment the room for labels
)
# this plots res.hcpc$call$X$clust
fviz_cluster(res.hcpc,
             cex = 0.8,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)
## Plot the most significant ones in each cluster (rank by test value only)
pull_cor <- function(df, name_list) {
    # dot to space, how ever there might be dots originally which will be replaced
    # wrongly. For example 1,7,7-trimethylbicyclo[2.2.1]hept-2-yl 4-hydroxy-3...
    # We could add those back but in general using names directly is not the best way.
    list_1 <- lapply(name_list,function(x) gsub("\\.", " ", x))
    # three space back to dot
    list_2 <- as.character(lapply(list_1,function(x) gsub("   ", "...", x)))
    # two space back to dot
    list_3 <- as.character(lapply(list_2,function(x) gsub("  ", ". ", x)))
    # get the coordinates
    df_cor <- data.frame(df[list_3,])
}
data_var <- res.pca$var$coord # the same as res.hcpc$call$t$res$var$coord
a <- data.frame(res.hcpc$desc.var$quanti[[1]])
b <- row.names(head(a[a$v.test>0,],10))
c <- data_var[b,]
data_var_1 <- data.frame(c,row.names = b)
data_var_1$cluster = paste0('cluster',1)

for (i in 2:length(res.hcpc$desc.var$quanti)) {
    print(i)
    a <- data.frame(res.hcpc$desc.var$quanti[[i]])
    b <- row.names(head(a[a$v.test>0,],5))
    c <- pull_cor(data_var,b)
    c$cluster <- paste0('cluster',i)
    data_var_1 <- rbind(data_var_1,c)
}

data_var_1$cluster <- as.factor(data_var_1$cluster)

ggplot(data_var_1) +
    geom_point(aes(x=Dim.1,y=Dim.2, color = cluster)) + 
    scale_color_jco() +
    geom_text_repel(aes(Dim.1,Dim.2,label = rownames(data_var_1))) +
    theme_minimal()

plot_NMDS <- function(data){
    #data_N <- data %>% select(-1) # select is also a function of MASS so be careful.
    N_matrix <- t(as.matrix(data))
    ##### NMDS
    # Collaps into two dimension using NMDS argorithm
    N_NMDS=metaMDS(N_matrix, # Our distance matrix
                   k=2) # The number of reduced dimensions
    # stress = 0.06 for negative Ions. <0.1 is considered great.
    # stress = 0.13 for Positive Ions. 
    
    stressplot(N_NMDS)
    # stress plot looks great too.
    # Now Plot 
    points <- as.data.frame(N_NMDS$points)
    colnames(points) <- c('NMDS1','NMDS2')
    g <- ggplot(points) +
        geom_point(aes(NMDS1,NMDS2),color = 'red') + 
        geom_text_repel(aes(NMDS1,NMDS2,label = rownames(points)))
    return(g)
}

plot_NMDS(data_F_sum)