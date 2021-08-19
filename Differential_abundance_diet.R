#differential abundance analyses
#example code shown for diet differential abundance analyses.

setwd("<>") #set working directory

#packages
library("DESeq2")
library("ggplot2")
library("dplyr")
library("phyloseq")
library("ashr") #to shrink log fold changes

#files
ps4<-readRDS("psWCE_phylo_18Oct20.rds")#contains 16s wild, captive, and enviro samples
psCacE5<-readRDS("diet1percent_rar_29May20.rds")  #diet data

#use unrarefied data
#remove already previously identified problem samples
Samples_toRemove <- c("S129_1138w", "S8_954w", "S32_992w", "S79_1051w", "S5_947w",
				"S131_1162w", "S133_1168w", "S48_1012w", "SW1302-D0")
ps.remove<-subset_samples(ps4, !(code %in% Samples_toRemove))
psWild<-subset_samples(ps.remove, Sample_type=="wild")
psWild20<-prune_taxa(taxa_sums(psWild) > 20, psWild) #keep taxa with > 20 reads

##alternative option to examine enriched families instead of ASVs
#psWild20<-tax_glom(psWild20, "Family")

#add column to sample data for each animal cac: y/n, creo: Y/N, cup: Y/N
#Y if at least 10% of diet
psCacE5<-subset_samples(psCacE5, code!="S79_1051w")
psF<-tax_glom(psCacE5, "Family")
merg = merge_samples(psF, "Site_Sp")
psFr = transform_sample_counts(merg, function(x) x/sum(x))

#Cactaceae
cacP<-subset_taxa(psFr, Family=="Cactaceae")
sumC<-data.frame(c=sample_sums(cacP), d=rep("cactus",25))
cacpops<-rownames(sumC[which(sumC$c>0.1),])

#Cupressaceae
cupP<-subset_taxa(psFr, Family=="Cupressaceae")
plot_bar(cupP)
sumCu<-data.frame(c=sample_sums(cupP), d=rep("Cup",25))
cuppops<-rownames(sumCu[which(sumCu$c>.1),])

#Zygophyllaceae
zP<-subset_taxa(psFr, Family=="Zygophyllaceae")
plot_bar(zP)
sumZ<-data.frame(c=sample_sums(zP), d=rep("zyg",25))
zpops<-rownames(sumZ[which(sumZ$c>.1),])

site_sp<-sample_data(psWild20)$Site_Sp

cacEat<-c(ifelse(site_sp %in% c(cacpops ), "Y", "N"))
cupEat<-c(ifelse(site_sp %in% c(cuppops ), "Y", "N"))
zygEat<-c(ifelse(site_sp %in% c(zpops ), "Y", "N"))

#add these as columns to the sample data
sample_data(psWild20)$cacEat<-as.factor(cacEat)
sample_data(psWild20)$cupEat<-as.factor(cupEat)
sample_data(psWild20)$zygEat<-as.factor(zygEat)

#use deseq2 to see if each diet (Y) is associated with increased abundance of
#taxa, compared to non-diet (N) feeders.
#levels (N/Y) ordered alphabetical such that N, null is first
#thus log2FoldChange > 0 means enriched in the Y plant feeders
#use shrunken LFC, repeat for each diet plant type
ps_dds<-phyloseq_to_deseq2(psWild20,design = ~ zygEat)  #cacEat, cupEat, zygEat

#DESeq cannot compute log geometric means when every gene contains zeros
#work around in https://github.com/joey711/phyloseq/issues/387
#calculate size factors using a zero-tolerant variant of geometric mean
#included in package vignette
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}
geoMeans = apply(counts(ps_dds), 1, gm_mean)
ps_dds = estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds = DESeq(ps_dds, fitType="local")

resAsh <- lfcShrink(ps_dds, coef="zygEat_Y_vs_N", type="ashr") #Change coef!!
summary(resAsh, alpha=0.05)
dfAsh<-as.data.frame(resAsh)

res<-dfAsh
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psWild20)[rownames(sigtab), ], "matrix"))
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ] #just increased

#repeat for each plant and save output:posigtab_cup, posigtab_creo, posigtab_cac
posigtab_creo<-posigtab
saveRDS(posigtab_creo, "pos_deseq_creo_20Oct20.R")

#For plots and downstream analyses
#merge into a single table, add diet type column first
posigtab_cac$diet<-as.factor(rep("cac", length(posigtab_cac$Genus)))
posigtab_cup$diet<-as.factor(rep("cup", length(posigtab_cup$Genus)))
posigtab_creo$diet<-as.factor(rep("zyg", length(posigtab_creo$Genus)))

posig_diets<-rbind(posigtab_cac,posigtab_cup, posigtab_creo)

####################################################
#Do these diet associated taxa decrease in captivity
#look at LFC of these taxa between wild and captive animals
#as above, but keeping both wild and captive sample_types
ps.r20<-prune_taxa(taxa_sums(ps.remove) > 20, ps.remove)

ps_dds<-phyloseq_to_deseq2(ps.r20,design = ~ Sample_type)  #

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}
geoMeans = apply(counts(ps_dds), 1, gm_mean)
ps_dds = estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds = DESeq(ps_dds, fitType="local")

res = results(ps_dds, cooksCutoff = FALSE)
sum(res$padj < 0.05, na.rm=TRUE)
res05 <- results(ps_dds, alpha=0.05) #set adj p to 0.05, instead of 0.1 default

resAsh <- lfcShrink(ps_dds, coef="Sample_type_wild_vs_captive", type="ashr")
summary(resAsh, alpha=0.05)
dfAsh<-as.data.frame(resAsh)

#the row names in posig_diets are the diet associated ASVs
dfAA<-merge(dfAsh, posig_diets, by=0)  #merge keeps LFC data for ASVs in both

#test whether diet associasted ASVs decrease in captivity
hist(dfAA$log2FoldChange.x)
mean(dfAA$log2FoldChange.x)
t.test(dfAA$log2FoldChange.x)
#################################################
