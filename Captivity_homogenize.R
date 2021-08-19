#does captivity homogenize the microbiome?

setwd("<>")#set file location

library(phyloseq)
library(vegan)

#files
ps.rout<-readRDS( "ps.rout_phylo_18Oct2020.rds")

psW<-subset_samples(ps.rout, Sample_type=="wild")
psC<-subset_samples(ps.rout, Sample_type=="captive")

#subset to animals with paired w/c data
pops_keep<-intersect(sample_data(psW)$Woodrat.ID, sample_data(psC)$Woodrat.ID)
ps.paired<-subset_samples(ps.rout, (Woodrat.ID %in% pops_keep)) 

bray_WC<-phyloseq::distance(ps.paired, "bray")#
modb <- betadisper(bray_WC, sample_data(ps.paired)$Sample_type)
modb
permutest(modb, pairwise = TRUE, permutations = 9999)

jac_WC<-phyloseq::distance(ps.paired, "jaccard")#
modj <- betadisper(jac_WC, sample_data(ps.paired)$Sample_type)
modj
permutest(modj, pairwise = TRUE, permutations = 9999)

un_WC<-phyloseq::distance(ps.paired, "unifrac")#
modu <- betadisper(un_WC, sample_data(ps.paired)$Sample_type)
modu
permutest(modu, pairwise = TRUE, permutations = 9999)

wu_WC<-phyloseq::distance(ps.paired, "wunifrac")#
modw <- betadisper(wu_WC, sample_data(ps.paired)$Sample_type)
modw
permutest(modw, pairwise = TRUE, permutations = 9999)

