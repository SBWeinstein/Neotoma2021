#differential abundance analyses, comparing wild and captive animals

setwd("<>") #set working directory

#packages
library("DESeq2")
library("ggplot2")
library("dplyr")
library("phyloseq")
library("ashr") #to shrink log fold changes
library("RColorBrewer")

#files
ps4<-readRDS("psWCE_phylo_18Oct20.rds")#contains wild, captive, and enviro samples

#use unrarefied wild and captive data
ps6<-subset_samples(ps4, Sample_type=="captive" |Sample_type=="wild")
#remove already previously identified problem samples
Samples_toRemove <- c("S129_1138w", "S8_954w", "S32_992w", "S79_1051w", "S5_947w",
				"S131_1162w", "S133_1168w", "S48_1012w", "SW1302-D0")
ps.remove<-subset_samples(ps6, !(code %in% Samples_toRemove))
#################################################
pops<-unique(sample_data(ps.remove)$Site_Sp) #get 25 population names
notpops<-c("WR_lep","JR_mac")  #remove pops that are wild only: WR_lep, JR_mac
popsWC<- setdiff(pops, notpops)
pops<-popsWC  #use these for loop

#run dif abundance for each pop, outputs list of dataframes, one for each pop
listofdfs <- list()

for(j in 1:length(pops)){
	pspop<-subset_samples(ps.remove, Site_Sp==pops[j]) #subset to a single population
	pspop<-prune_taxa(taxa_sums(pspop) > 0, pspop) #pop1 688 taxa
	pspop20<-prune_taxa(taxa_sums(pspop) > 20, pspop) #keep taxa with > 20 reads, 509 taxa in pop1

	pop_dds<-phyloseq_to_deseq2(pspop20,design = ~Sample_type  )  #comparing wild v captive

	gm_mean = function(x, na.rm=TRUE){
 		 exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
		}
	geoMeans = apply(counts(pop_dds), 1, gm_mean)
	pop_dds = estimateSizeFactors(pop_dds, geoMeans = geoMeans)
	pop_dds = DESeq(pop_dds, fitType="local")

	resAsh <- lfcShrink(pop_dds, coef="Sample_type_wild_vs_captive", type="ashr")
	#res = results(pop_dds, cooksCutoff = FALSE)
	alpha = 0.05
	sigtab = resAsh[which(resAsh$padj < alpha), ]
	sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pspop20)[rownames(sigtab), ], "matrix"))

	sigtab$pop<-as.factor(rep(pops[j], length(sigtab$Genus)))
	sigtab$taxa<-rownames(sigtab)
	rownames(sigtab) <- c()
	listofdfs[[j]] <- sigtab
	}

allpops2<-bind_rows(listofdfs)  #bind all of the the dataframes into one based on columns

#the long taxa names are hard to look at. give them short codes
sigtaxa<-unique(allpops2$taxa)
prefix<-"ASV"
nums<-seq(1:length(sigtaxa))
codes<- paste (prefix, nums, sep= "_")
taxa_codes<-data.frame(ASV=sigtaxa, codes=codes)
allpops3<-merge(allpops2, taxa_codes, by.x= "taxa", by.y="ASV")

#add column for direction of change
allpops3$dir<-as.factor(ifelse(allpops3$log2FoldChange>0, "P", "N"))

#for ggplot, create long form data frame with pop, taxa codes, and log2foldchange
dlong<-data.frame(pop=allpops3$pop, taxa=allpops3$taxa,
	                L2Ch=allpops3$log2FoldChange, phylum=allpops3$Phylum )

#how many asvs change per pop?
sig_cts1<-table(allpops2$pop)
mean(sig_cts1)
sd(sig_cts1)

#how many up v down? (up is enriched in wild)
#levels (C/W) ordered alphabetical such that C, null is first
#thus log2FoldChange" > 0 means enriched in the W
sig_cts2<-table(allpops3$pop, allpops3$dir)
sig_cts3<-as.data.frame.matrix(sig_cts2)  #table to data frame, finicky

mean(sig_cts3$N)
sd(sig_cts3$N)

mean(sig_cts3$P) #enriched in wild (lower in captivity)
sd(sig_cts3$P)

sig_cts4<-table(allpops3$codes)
mean(sig_cts4)
sd(sig_cts4)
hist(sig_cts4)
table(sig_cts4)#369 occur in only one pop
dim(sig_cts4)#475 ASVs
369/475 # 0.6849817 of ASVs detected in one pop

#what are most taxa
fams<-table(allpops3$Family)
229/475 #.48, 229 of 475 taxa in Muribaculaceae
149/475 #.31, 149 of 475 taxa in Lachnospiraceae

#how many taxa are M or L in dataset?
ps.remove  #W and C- 2754 taxa
psM<-subset_taxa(ps.remove, Family== "Muribaculaceae")#1187,1187/2754  = 0.43
psL<-subset_taxa(ps.remove, Family== "Lachnospiraceae")#707, 707/2754  = 0.26

#####
#make plots 6A,B for appendix
# Phylum order
x = tapply(allpops3$log2FoldChange, allpops3$Phylum, function(x) max(x))
x = sort(x, TRUE)
allpops3$Phylum = factor(as.character(allpops3$Phylum), levels=names(x))
# Family order
x = tapply(allpops3$log2FoldChange, allpops3$Family, function(x) max(x))
x = sort(x, TRUE)
allpops3$Family = factor(as.character(allpops3$Family), levels=names(x))

#add pop short codes
cds2<-read.csv("Site_Sp_codes.csv")
cds2<-data.frame(code = cds2$code2, Site_Sp=cds2$Site_Sp)
allpops4<-merge(allpops3, cds2, by.x=13, by.y=2)
allpops4N<-allpops4[which(allpops4$dir=="N"),]
allpops4P<-allpops4[which(allpops4$dir=="P"),]

theme_set(theme_bw())
#for A and B use allpops4N or allpops4P
ggplot(allpops4P, aes(x=Family, y=log2FoldChange, color=Family)) +
	geom_jitter(width=.2) +
	facet_grid(code~., space="free")+
	scale_y_continuous(breaks=c(0, 10,20,30),labels=c(0,10,20, 30))+  #For P (more abundant in wild)
	#scale_y_continuous(breaks=c(5, -10, -25),labels=c(5, -10, -25))+  #For N (more abundant in cap)
  theme(axis.text.x = element_text(angle = -75, hjust = 0, vjust=0.5), legend.position = "none")

##############################
#across all pops, are families (or genera) increasing/decreasing?
#also remove pops that are wild only
ps.paired<-subset_samples(ps.remove, !(Site_Sp %in% notpops))
psp1<-prune_taxa(taxa_sums(ps.paired) > 0, ps.paired)

pspF<-tax_glom(psp1,"Genus" ) #for Genus: 95 taxa, Repeat for "Family": 47 taxa

pF_dds<-phyloseq_to_deseq2(pspF,design = ~Sample_type  )  #comparing wild v captive
	gm_mean = function(x, na.rm=TRUE){
 		 exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
		}
	geoMeans = apply(counts(pF_dds), 1, gm_mean)
	pF_dds = estimateSizeFactors(pF_dds, geoMeans = geoMeans)
	pF_dds = DESeq(pF_dds, fitType="local")
	#shrink LFC
	res <- lfcShrink(pF_dds, coef="Sample_type_wild_vs_captive", type="ashr")
	alpha = 0.05
	sigtab = res[which(res$padj < alpha), ]
	sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pspF)[rownames(sigtab), ], "matrix"))

#in sigtab, rownames are sequences
df2<-data.frame(Family=sigtab$Family,Genus=sigtab$Genus, log2FoldChange=sigtab$log2FoldChange)
df2$dir<-as.factor(ifelse(df2$log2FoldChange>0, "P", "N"))
df2<-df2[order(df2$dir),]
#for ggplot, create long form data frame with pop, taxa codes, and log2foldchange
dlong<-data.frame( taxa=df2$Family, L2Ch=df2$log2FoldChange, pop=rep(1, length(df2$Family)))

#############################
#make heat map for figure 5D
heatmap.plot <- ggplot(data = dlong, aes(x = reorder(taxa,-L2Ch), y=pop)) +
  geom_tile(aes(fill = L2Ch)) +
 scale_fill_distiller(palette = "BrBG",direction=1)  +
	theme_classic()+
  theme(axis.text.y = element_blank(),
		axis.title = element_blank(),
		axis.ticks.y = element_blank(),
		 axis.text.x=element_text(angle=60, hjust=1))
