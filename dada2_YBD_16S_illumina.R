library(dada2)
library(ggplot2)
library(phyloseq)
library(phangorn)
library(DECIPHER)

path <- '~/Data/YBD/16S_illumina/fastq/'
files = list.files(path, pattern = '.fastq')

# we also need the sample names
sample_names <- sapply(strsplit(files,"\\."),
                       `[`,  # backtick (or back quote, `), [ for exact match while [[ can extend
                       1) # indexing by character, extracts the first element of a subset

trim_files = list.files('~/Data/YBD/16S_illumina/trimmed/', pattern = '.fastq', full.names = TRUE)

for (f in files){
    print(f) # Reserved words are not valid identifiers. file is not one of them. but somehow it is taken before using.
    plotQualityProfile(paste0(path,f))
}

# need to filter the long&crappy reads (using python), later
out <- filterAndTrim(files, trim_files, truncLen=290, maxN=0,
                     maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE,
                     multithread=TRUE)
# 
errors <- learnErrors(trim_files, multithread=TRUE)
# time demanding. But why only 12 samples

derep <- derepFastq(trim_files, verbose=TRUE)

# name the derep-class objects by the sample names
names(derep) <- sample_names

dada_1 <- dada(derep, err=errors, multithread=TRUE)
seq_table <- makeSequenceTable(dada_1)
dim(seq_table)

# inspect distribution of sequence lengths
table(nchar(getSequences(seq_table)))

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
                       '~/Data/MiSeq_SOP/silva_nr_v128_train_set.fa.gz',
                       multithread=TRUE)

# time consuming

seq_table <- readRDS('seq_table.rds')
taxa <- addSpecies(taxa, '~/Data/MiSeq_SOP/silva_species_assignment_v128.fa.gz')
saveRDS(taxa, 'taxa.rds')
saveRDS(seq_table, 'seq_table.rds')
saveRDS(seq_table_nochim, 'seq_table_nochim.rds')
taxa <- readRDS('taxa.rds')
taxa_print <- taxa  # removing sequence rownames for display only
rownames(taxa_print) <- NULL
head(taxa_print)

sequences <- getSequences(seq_table)
names(sequences) <- sequences  # this propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
saveRDS(phang_align,'phang_align.rds')
dm <- dist.ml(phang_align)
treeNJ <- NJ(dm)  # note, tip order != sequence order
# takes a long time! Also I don't need this thing.
saveRDS(treeNJ, 'treeNJ.rds')
treeNJ <- readRDS('treeNJ.rds')
phang_align <- readRDS('phang_align.rds')
# pml computes the likelihood of a phylogenetic tree given a sequence alignment and a model. optim.pml optimizes the different model parameters.
fit = pml(treeNJ, data=phang_align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
saveRDS(fitGTR, 'fitGTR.rds')
fitGTR <- readRDS('fitGTR.rds')
fitGTR <- optim.pml(fitGTR, model='GTR', optInv=TRUE, optGamma=TRUE,
                    rearrangement = 'stochastic',
                    control = pml.control(trace = 0))
#Error: vector memory exhausted (limit reached?)
# First, check whether you are using 64-bit R. 
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
