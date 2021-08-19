#Run in 64 bit,  R version (3.6.3)
#code based on:
#Dada2: https://benjjneb.github.io/dada2/tutorial.html

###Load packages
.cran_packages <- c("ggplot2", "gridExtra", "knitr")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn", "BiocStyle")
# Load packages into session
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
library("ShortRead")
##################################################
#set working directory
setwd("<directory_location>") #input path here
#directory with fastq files.
path <- "./<Directory_with_fastq.gz_files_from_1_run>"

#tell R where to find cutadapt
cutadapt <- "<cutadapt_location>" #cutadapt path
system2(cutadapt, args = "--version") # Run shell commands from R

###############from FastQ files to (merged) OTU table...#####################
#use Cutadapt to first remove primers and filter sequences
#generate matched lists of the forward and reverse read files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

FWD <- "GGGCAATCCTGAGCCAA"  ##  forward primer sequence, trnl-g
REV <- "CCATTGAGTCTCTGCACCTATC"  ## Reverse primer sequence, trnl-h

#verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna),
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))
  }
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#Remove sequences with ambiguous bases.
#Put N-filtered files in filtN/ subdirectory
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)

##count the number of times the primers appear in the forward and reverse reads
primerHits <- function(primer, fn) {
  #Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
  }
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                              "--discard-untrimmed", "--minimum-length", 8,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                             fnFs.filtN[i], fnRs.filtN[i]))
}

#sanity check,  count the presence of primers in a cutadapt-ed sample (#2)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[2]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[2]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[2]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[2]]))

#primer-free sequence files are now ready to be analyzed in DADA2 pipeline
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#inspect read quality
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#filter and trim sequences
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
      truncQ = 2, minLen = 10, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
head(out)

##estimate error rates
#pcr blank in run464 had zero sequences after filtering
#remove samples with no remaining sequences before proceeding
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#infer sequence variants using DADA2
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
#dereplicate data
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)
#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

seqtab4 <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab4)))

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)
track

#save sequence table from each run
saveRDS(seqtab4, "<file_name>.rds")
#####################################################
#Repeat cutadapt, and sample inference for each run seperatly
#################################################
#merge  sequence tables ("464","435","MarchAB")
st.all <- mergeSequenceTables(R435, R464, MarchAB)
#remove chimeras from merged sequence table
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=FALSE)
write.csv(seqtab.nochim, "seqtab_trnL_nochim_5April20.csv")

#Use plant sequence table CSV for taxonomy assignment using python script
