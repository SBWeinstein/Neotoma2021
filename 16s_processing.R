#Initial processing of 16s amplicon sequences
#starts with fastq files, produces single phyloseq object for the project
#Run on 64 bit, Redhat linux system, Rstudio 3.6.1

#code based on:
#Callahan et al 2016 Pipeline: https://f1000research.com/articles/5-1492
#Dada2: https://benjjneb.github.io/dada2/tutorial.html
#combining sequencing runs: https://benjjneb.github.io/dada2/bigdata.html
# using cutadapt: https://benjjneb.github.io/dada2/ITS_workflow.html
#Phyloseq:https://joey711.github.io/phyloseq/

###Load packages
.cran_packages <- c("ggplot2", "gridExtra", "knitr", "ShortRead", "ips")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn", "BiocStyle")
# Load packages into session
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

###############################################
#set working directory
setwd("<directory_location>") #input working directory here
#location of sequence files for each illumina run
pathA <- "./<Directory_with_fastq.gz_files_from_1_run>" #input file name here

#tell R where to find cutadapt
cutadapt <- "<cutadapt_location>" #input path here
system2(cutadapt, args = "--version") # Run shell commands from R

#Other files:
#sample metadata
samdf <- read.csv("metadata_16s_phylosym_2April20.csv", header=TRUE)
##############################################
#from FastQ files to (merged) OTU table
#16s  sequences are from 4 runs
#run cutadapt and sample inference script independently for
#each run to account for diffferent error rates.
#Then merge into a full study table before chimera removal/taxonomy assignments

###use Cutadapt to first remove primers and filter sequences
#make matched lists of forward and reverse read files, parse out sample names
fnFs <- sort(list.files(pathA, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(pathA, pattern="_R2_001.fastq", full.names = TRUE))

FWD <- "GTGYCAGCMGCCGCGGTAA"  ##  forward primer sequence, 515F, updated
REV <- "GGACTACNVGGGTWTCTAAT"  ## Reverse primer sequence, 806R, updated

#verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna),
                 Reverse = reverse(dna),
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  #Convert back to character vector
  }
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#Remove sequences with ambiguous bases.
#Put N-filtered files in filtN/ subdirectory
fnFs.filtN <- file.path(pathA, "filtN", basename(fnFs))
fnRs.filtN <- file.path(pathA, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#count the number of times primers appear in the forward and reverse read,
#considering all possible primer orientations.
primerHits <- function(primer, fn) {
  #Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#use cutadapt to trim primers
path.cut <- file.path(pathA, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
#Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
#Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                              "--discard-untrimmed", "--minimum-length", 8,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
                           }

#sanity check,  count the presence of primers in a cutadapt-ed sample
sq=2 #pick a sample
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[sq]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[sq]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[sq]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[sq]]))

#primer-free sequence files are ready for DADA2 pipeline
#Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#inspect read quality
plotQualityProfile(cutFs[8:9])
plotQualityProfile(cutRs[8:9])

####################Read filtering and Trimming
#Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#filter and trim sequences
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out) #look at first 5

#estimate error rates
errF <- learnErrors(filtFs, multithread =FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Dereplication
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

#Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#dereplicate data
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab16s_Run2 <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab16s_Run2)))

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)

#save sequence table from each run
saveRDS(seqtab16s_Run2, "seqtab_16s_Run2.rds")
saveRDS(track, "read_tracker_Run2.rds")

###############################################################################
#Repeat cutadapt, and sample inference for each sequencing run
#merge 4 sequence tables into one (Run1, Run2, Run3, RunAB)
seqtab16s_merged <- mergeSequenceTables(seqtab16s_RunMAB, seqtab16s_Run1, seqtab16s_Run2, seqtab16s_Run3)
#remove chimeras from merged sequence table
seqtab.nochim <- removeBimeraDenovo(seqtab16s_merged, method="consensus", multithread=TRUE)
#fraction of merged sequences retained after chimera removal (~97%)
sum(seqtab.nochim)/sum(seqtab16s_merged)
#save merged, no chimera table for downstream analyses.
saveRDS(seqtab.nochim, "seqtab_merged_nochim_1April20.rds")

###############################################################################
#assign taxonomy, as of 1April2020, current silva release is v138
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v138_train_set.fa.gz", multithread=FALSE)
saveRDS(taxa, "taxa_merged_nochim_1April20.rds")

#construct phylogenetic tree
#export aligned sequences and use FastTree
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
#Multiple Sequence Alignment via DECIPHER
mult <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
#write out alignment file in FASTA format for FastTree in Terminal
#use writeFASTA (shortreads) and not write.FASTA (ape)
writeFasta(mult, "aligned_seqs.fasta")
#use FastTree to produce "MLtree2.nwk"
#code (RUN IN FastTree. NOT IN R): FastTree -gtr -nt < msa2.fasta > MLtree2.nwk

##############################################################################
#Make phyloseq object
#read in csv of sample data, already in same order as sequence table data
samdf <- read.csv("metadata_16s_phylosym_2April20.csv", header=TRUE)
#add a new column with the sample names used in the sequence table,
#set first colum to be row names
samdf1<-cbind(rownames(seqtab.nochim), samdf)
samdf2 <- data.frame(samdf1, row.names = 1)
#upload ML tree from FastTree
tree<-read.tree("MLtree2.nwk")

#make phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf2),
               tax_table(taxa), phy_tree(tree))

#this contains some samples for other projects, keep PhyloProj == Y or q
ps1<-subset_samples(ps, PhyloProj=="Y"|PhyloProj=="q")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)

#Exploration of environmental sources, blanks, etc not included here
#remove blanks
ps3<-subset_samples(ps2, Sample_type=="captive" |
                         Sample_type=="wild"|
                         Sample_type=="environ_control")
ps3<-prune_taxa(taxa_sums(ps3) > 0, ps3)

#look at prevalence and abundance of taxa to make a filtering decision
#Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps3),
                 MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
#Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps3),
                      tax_table(ps3))

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps3),color=Phylum)) +
  geom_hline(yintercept = 0.03, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") +
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#Define prevalence threshold as 2% of total samples
prevalenceThreshold = 0.02 * nsamples(ps3)

# Execute prevalence filter, using prune_taxa()
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps4 = prune_taxa(keepTaxa, ps3)

#2% 2754 taxa in 5489658
sum(sample_sums(ps4))

#contains wild, captive, and enviro samples
saveRDS(ps4, "psWCE_phylo_18Oct20.rds")

#remove environmental controls
ps6<-subset_samples(ps4, Sample_type=="captive" |Sample_type=="wild")

ps.rar = rarefy_even_depth(ps6, rngseed=1, sample.size=5300, replace=F)

#quality control, check for outliers and examine original DNA quality
#outliers to remove- (wunifrac,W&C) 1138 wild, 954 wild extreme outliers, poor dna quality;
	#1051w,947w poor DNA quality, unlike others from pop or species
	#992W- suspected wrong sample was extracted
	#1162, 1168, 1012-(wilds) poor dna quality, unusual taxa distribution in bar_plots
	#1302- poor amplification (already removed by rarefy)

Samples_toRemove <- c("S129_1138w", "S8_954w", "S32_992w", "S79_1051w", "S5_947w",
				"S131_1162w", "S133_1168w", "S48_1012w", "SW1302-D0")
ps.rout<-subset_samples(ps.rar, !(code %in% Samples_toRemove))

saveRDS(ps.rout, "ps.rout_phylo_18Oct2020.rds")
