# https://sgbc.github.io/course/tutorials/docs/16S/#metabarcoding
# Install and Load Packages
# R version 4.0.2
# source('https://bioconductor.org/biocLite.R') # obsolete
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
# sudo apt-get install libbz2-dev
BiocManager::install("Rhtslib")
BiocManager::install("Rsamtools")
BiocManager::install('GenomicAlignments')
BiocManager::install('codetools')
BiocManager::install('ShortRead')
install.packages('Rcpp')
BiocManager::install('dada2')
install.packages('lattice')
BiocManager::install('phyloseq')
install.packages('~/Downloads/sf_0.9-6.tgz', repos = NULL, type="binary")
install.packages('~/Downloads/rgdal_1.5-18.tgz', repos = NULL, type="binary")
BiocManager::install('DECIPHER')
install.packages('ggplot2')
install.packages('phangorn')
BiocManager::install("Rhdf5lib")

library(dada2)
library(ggplot2)
library(phyloseq)
library(phangorn)
library(DECIPHER)
packageVersion('dada2')

path <- '~/Data/MiSeq_SOP'
list.files(path)

# Filtering and Trimming
# First we create two lists with the sorted name of the reads: one for the forward reads, one for the reverse reads

raw_forward <- sort(list.files(path, pattern="_R1_001.fastq",
                               full.names=TRUE))

raw_reverse <- sort(list.files(path, pattern="_R2_001.fastq",
                               full.names=TRUE))

# we also need the sample names
sample_names <- sapply(strsplit(basename(raw_forward), "_"),
                       `[`,  # extracts the first element of a subset
                       1)

# then we visualise the quality of our reads
plotQualityProfile(raw_forward[1:2]) # the first two
plotQualityProfile(raw_reverse[1:2]) 
plotQualityProfile(c(raw_forward[1],raw_reverse[1]))

# place filtered files in filtered/ subdirectory
filtered_path <- file.path(path, "filtered")

filtered_forward <- file.path(filtered_path,
                              paste0(sample_names, "_R1_trimmed.fastq.gz"))

filtered_reverse <- file.path(filtered_path,
                              paste0(sample_names, "_R2_trimmed.fastq.gz"))

out <- filterAndTrim(raw_forward, filtered_forward, raw_reverse,
                     filtered_reverse, truncLen=c(240,160), maxN=0,
                     maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE,
                     multithread=TRUE)
head(out)

errors_forward <- learnErrors(filtered_forward, multithread=TRUE)
errors_reverse <- learnErrors(filtered_reverse, multithread=TRUE)

derep_forward <- derepFastq(filtered_forward, verbose=TRUE)
derep_reverse <- derepFastq(filtered_reverse, verbose=TRUE)
# name the derep-class objects by the sample names
names(derep_forward) <- sample_names
names(derep_reverse) <- sample_names

dada_forward <- dada(derep_forward, err=errors_forward, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=errors_reverse, multithread=TRUE)

# inspect the dada-class object
dada_forward[[1]]

merged_reads <- mergePairs(dada_forward, derep_forward, dada_reverse,
                           derep_reverse, verbose=TRUE)

# inspect the merger data.frame from the first sample
head(merged_reads[[1]])

seq_table <- makeSequenceTable(merged_reads)
dim(seq_table)

# inspect distribution of sequence lengths
table(nchar(getSequences(seq_table)))

seq_table_nochim <- removeBimeraDenovo(seq_table, method='consensus',
                                       multithread=TRUE, verbose=TRUE)
dim(seq_table_nochim)

# which percentage of our reads did we keep?
sum(seq_table_nochim) / sum(seq_table)

get_n <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dada_forward, get_n), sapply(merged_reads, get_n),
               rowSums(seq_table), rowSums(seq_table_nochim))

colnames(track) <- c('input', 'filtered', 'denoised', 'merged', 'tabled',
                     'nonchim')
rownames(track) <- sample_names
head(track)

taxa <- assignTaxonomy(seq_table_nochim,
                       '~/Data/MiSeq_SOP/silva_nr_v128_train_set.fa.gz',
                       multithread=TRUE)
taxa <- addSpecies(taxa, '~/Data/MiSeq_SOP/silva_species_assignment_v128.fa.gz')

taxa_print <- taxa  # removing sequence rownames for display only
rownames(taxa_print) <- NULL
head(taxa_print)

sequences <- getSequences(seq_table)
names(sequences) <- sequences  # this propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
dm <- dist.ml(phang_align)
treeNJ <- NJ(dm)  # note, tip order != sequence order
fit = pml(treeNJ, data=phang_align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model='GTR', optInv=TRUE, optGamma=TRUE,
                    rearrangement = 'stochastic',
                    control = pml.control(trace = 0))
detach('package:phangorn', unload=TRUE)

sample_data <- read.table(
    'https://hadrieng.github.io/tutorials/data/16S_metadata.txt',
    header=TRUE, row.names="sample_name")

physeq <- phyloseq(otu_table(seq_table_nochim, taxa_are_rows=FALSE),
                   sample_data(sample_data),
                   tax_table(taxa),
                   phy_tree(fitGTR$tree))
# remove mock sample
physeq <- prune_samples(sample_names(physeq) != 'Mock', physeq)
physeq
plot_richness(physeq, x='day', measures=c('Shannon', 'Fisher'), color='when') +
    theme_minimal()

ord <- ordinate(physeq, 'MDS', 'euclidean')
plot_ordination(physeq, ord, type='samples', color='when',
                title='PCA of the samples from the MiSeq SOP') +
    theme_minimal()

ord <- ordinate(physeq, 'NMDS', 'bray')
plot_ordination(physeq, ord, type='samples', color='when',
                title='PCA of the samples from the MiSeq SOP') +
    theme_minimal()

top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]
physeq_top20 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
physeq_top20 <- prune_taxa(top20, physeq_top20)
plot_bar(physeq_top20, x='day', fill='Family') +
    facet_wrap(~when, scales='free_x') +
    theme_minimal()

bacteroidetes <- subset_taxa(physeq, Phylum %in% c('Bacteroidetes'))
plot_tree(bacteroidetes, ladderize='left', size='abundance',
          color='when', label.tips='Family')
