#test for phylogenetic signal in captivity induced changes
#in diversity and composition for host species
#packages phyloseq and phylosignal have name space conflicts
#simpler to run from a clean workspace

setwd("C:/Users/Sara//Dropbox/phylosymbiosis/Sequence_analyses_Mar2020/16s_analyses")
capdf2<-readRDS("capdf2_4Aug21.rds")
Neo_tree<-read.tree(file="C:/Users/Sara//Dropbox/phylosymbiosis/neotoma/Neotoma_phylosymtree.nwk")

library(phangorn)#for reading in tree
library(phylobase)  #for function phylo4d
library(phylosignal)  #for Morans I

#get average values for each species for Per.Change.Obs
#using dataframe from Captivity_impacts.R script
pco<-data.frame(tapply(capdf2$Per.Change.Obs, capdf2$Species, mean))
names(pco)<-c("od")  #rename for observed diff
br<-tapply(capdf2$bray, capdf2$Species, mean)
ja<-tapply(capdf2$Jac, capdf2$Species, mean)
wu<-tapply(capdf2$Wunif, capdf2$Species, mean)
u<-tapply(capdf2$Unif, capdf2$Species, mean)

mns<-data.frame(pco, br,ja,wu,u)

#replace tip labels to just species
Neo_tree$tip.label<-c("N. albigula", "N. stephensi", "N. macrotis", 
                      "N. cinerea", "N. bryanti", "N. lepida", "N. devia" )
Neo_tree$node.label<-c(letters[1:6])#phylo4d needs unique node names
p4d <- phylo4d(Neo_tree, mns) #merge data and tree
#visualize patterns
barplot.phylo4d(p4d, tree.type = "phylo", tree.ladderize = FALSE)
sig<-phyloSignal(p4d = p4d, method = "I", reps = 999)
sigs<-c(sig$pvalue)

p.adjust(rev(sigs$I), method="fdr")#bonferroni or fdr

#not significant
#also tested at population level, also not significant, not shown


