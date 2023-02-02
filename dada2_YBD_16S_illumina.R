library(dada2)
library(ggplot2)
library(phyloseq)
library(phangorn)
library(DECIPHER)
library(dplyr)

path <- '~/Data/YBD/16S_illumina/fastq/'
files = list.files(path, pattern = '.fastq', full.names = TRUE)
# plotQualityProfile(files)
# we also need the sample names
files_list = list.files(path, pattern = '.fastq')
file_names <- sapply(strsplit(files_list,"\\."),
                       `[`,  # backtick (or back quote, `), [ for exact match while [[ can extend
                       1) # indexing by character, extracts the first element of a subset
sample_data <- read.csv('~/Data/YBD/YBD_metadata.csv',
                        header=TRUE)
sample_names <- sample_data[match(file_names,sample_data$ID),]$Sample
rownames(sample_data) <- sample_data$Sample
filtered_path <- '~/Data/YBD/16S_illumina/trimmed'
trim_files <- file.path(filtered_path,
                              paste0('trim_', sample_names, ".fastq.gz"))



# Primer used: 515F-806R, 292bp, V4 of 16S
# need to filter the long&crappy reads (using python), later
out <- filterAndTrim(files, trim_files, maxLen = 300, maxN=0,
                     maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE,
                     multithread=TRUE, trimLeft = 23, trimRight = 20)
# 
errors <- learnErrors(trim_files, multithread=TRUE)
saveRDS(errors,'~/Data/YBD/16S_illumina/dada2_errors.rds')
# time demanding. But why only 12 samples

derep <- derepFastq(trim_files, verbose=TRUE)

# name the derep-class objects by the sample names
names(derep) <- sample_names

dada_1 <- dada(derep, err=errors, multithread=TRUE)
seq_table <- makeSequenceTable(dada_1)
dim(seq_table)

# inspect distribution of sequence lengths
table(nchar(getSequences(seq_table)))
#seqtab2 <- seq_table[,nchar(colnames(seq_table)) %in% 250:256]
seq_table_nochim <- removeBimeraDenovo(seq_table, method='consensus',
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
                       '~/Data/silva_nr_v128_train_set.fa.gz',
                       multithread=TRUE)

# time consuming

#seq_table <- readRDS('seq_table.rds')
taxa <- addSpecies(taxa, '~/Data/silva_species_assignment_v128.fa.gz')
saveRDS(taxa, '~/Data/YBD/16S_illumina/dada2_taxa.rds')
#saveRDS(seq_table, 'seq_table.rds')
#saveRDS(seq_table_nochim, 'seq_table_nochim.rds')
#taxa <- readRDS('taxa.rds')
taxa_print <- taxa  # removing sequence rownames for display only
rownames(taxa_print) <- NULL
head(taxa_print)

sequences <- getSequences(seq_table)
names(sequences) <- sequences  # this propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA) # time consuming
phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
saveRDS(phang_align,'~/Data/YBD/16S_illumina/dada2_phang_align.rds')
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
OTU4 <- OTU4 %>% select(OTU, Cow_Thi_2:Cow_Sav_5)
write.table(OTU4,'OTU4.csv',row.names = FALSE, sep = ',', quote = FALSE)

taxa_names(physeq1) <- paste0(">SV", seq(ntaxa(physeq1)))
OTU5 <- as(otu_table(physeq1), "matrix")
OTU6 <- as.data.frame(t(OTU5))
OTU6$Seq <- colnames(OTU1)
OTU6 <- OTU6 %>% select(Seq)
write.table(OTU6, file = 'OTU6.fa', quote = FALSE, sep = '\n', col.names = FALSE)


plot_richness(physeq, x='Animal', measures=c('Shannon', 'Fisher'), color='Habitat') +
    theme_minimal()

ord <- ordinate(physeq, 'MDS', 'euclidean')
plot_ordination(physeq, ord, lable='Sample', color='Animal', shape = 'Habitat', 
                title='PCA of 16S data') +
    theme_minimal()

ord <- ordinate(physeq, 'NMDS', 'bray')
plot_ordination(physeq, ord, type='samples', color='Species', shape = 'Habitat',
                title='PCA of the samples from the MiSeq SOP') +
    theme_minimal()

top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]
physeq_relative <- transform_sample_counts(physeq1, function(OTU) OTU/sum(OTU))
physeq_top20 <- prune_taxa(top20, physeq_relative)

plot_bar(physeq_top20, x="ID", fill='Family') +
    facet_wrap(~Habitat, scales='free_x') +
    theme_minimal()

plot_bar(physeq_top20, x="ID", fill='Family') +
    facet_wrap(~Animal, scales='free_x') +
    theme_minimal()

physeq1 <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE),
                    sample_data(sample_data1),
                    tax_table(taxa))
physeq_top20_count <- prune_taxa(top20, physeq1)
OTU1 <- as(otu_table(physeq_top20_count), "matrix")
OTU2 <- t(OTU1)
taxa_names(physeq_top20) <- paste0("SV", seq(ntaxa(physeq_top20)))
OTU3 <- as(otu_table(physeq_top20), "matrix")
OTU4 <- as.data.frame(t(OTU3))
OTU4$OTU <- rownames(OTU4)
OTU4 <- OTU4 %>% select(OTU, A:R)
write.table(OTU4,'OTU4.csv',row.names = FALSE, sep = ',', quote = FALSE)
taxa_names(physeq_top20) <- paste0(">SV", seq(ntaxa(physeq_top20)))
OTU5 <- as(otu_table(physeq_top20), "matrix")
OTU6 <- as.data.frame(t(OTU5))
OTU6$Seq <- colnames(OTU1)
OTU6 <- OTU6 %>% select(Seq)
write.table(OTU6, file = 'OTU6.fa', quote = FALSE, sep = '\n', col.names = FALSE)

colnames(OTU1)
for (i in 1:nrow(OTU1)){
    
}