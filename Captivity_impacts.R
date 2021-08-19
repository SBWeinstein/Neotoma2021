#impacts of captivity

library(ggplot2)
library(phyloseq)
library(vegan)
library(tidyr) #for gather()
library(emmeans) #posthoc model comparisons

setwd("<>") #set working directory

#files
ps.rout<-readRDS("ps.rout_phylo_18Oct2020.rds") #microbiome data
plant<-readRDS( "diet1percent_rar_29May20.rds") #diet data
site_codes<-read.csv("Site_Sp_codes.csv") #site short codes for plots

#do captive animals retain population signal?
psC<-subset_samples(ps.rout, Sample_type=="captive") #subset to captive
sampledfC <- data.frame(sample_data(psC))

psC_bray<-phyloseq::distance(psC, "bray")
psC_wu<-phyloseq::distance(psC, "wunifrac")
psC_jac<-phyloseq::distance(psC, "jaccard")
psC_u<-phyloseq::distance(psC, "unifrac")

# Adonis test
adonis(psC_bray ~ Site_Sp, data = sampledfC)  #sig
adonis(psC_jac ~ Site_Sp, data = sampledfC)  #sig
adonis(psC_wu ~ Site_Sp, data = sampledfC)  #sig
adonis(psC_u ~ Site_Sp, data = sampledfC)  #sig

#does community composition differ between wild and captive animals?
# make a data frame from the sample_data
df_rt <- data.frame(sample_data(ps.rout))

psWC_bray<-phyloseq::distance(ps.rout, "bray")
psWC_wu<-phyloseq::distance(ps.rout, "wunifrac")
psWC_jac<-phyloseq::distance(ps.rout, "jaccard")
psWC_u<-phyloseq::distance(ps.rout, "unifrac")

# Adonis test
adonis(psWC_bray ~ Sample_type, data = df_rt)  #sig
adonis(psWC_jac ~ Sample_type, data = df_rt)  #sig
adonis(psWC_wu ~ Sample_type, data = df_rt)  #sig
adonis(psWC_u ~ Sample_type, data = df_rt)  #sig

#to examine factors that influence alpha and beta diversity metrics
#use only animals with both wild and captive data
#create data frame with animal info, diversity metrics, etc
psW<-subset_samples(ps.rout, Sample_type=="wild")#use wild samples to start
df <- data.frame(estimate_richness(psW, measures = c("Observed", "Shannon") ),
	Woodrat.ID=sample_data(psW)$Woodrat.ID,
	Species=sample_data(psW)$Species,
	site=sample_data(psW)$site,
	Site_Sp=sample_data(psW)$Site_Sp,
	conc=sample_data(psW)$ng.ul,
	run=sample_data(psW)$Run )

##add columns for diet observed families, and diet shannon diversity
#taxa glom at family level
psF<-tax_glom(plant, "Family")
#observed, shannon richness of families per sample
df_plant<-data.frame(estimate_richness(psF, measures=c("Observed", "Shannon")),
				Woodrat.ID=sample_data(psF)$Woodrat.ID)

#merge plant and 16s wild info, keep rats with missing plant data
dfPB<-merge(df, df_plant, by= "Woodrat.ID", all.x=TRUE)
colnames(dfPB)[c(2,3,9,10)] <-c("Obs.16sW","Shan.16sW", "Obs.plant", "Shan.plant")

#get observed and shannon 16s info for captive rats ("psC")
dfC <- data.frame(estimate_richness(psC, measures = c("Observed", "Shannon") ),
	Woodrat.ID=sample_data(psC)$Woodrat.ID)

#merge by woodrat Id, returns dataframe of only individuals with wild and cap
dfPBC<-merge(dfPB, dfC, by= "Woodrat.ID")
colnames(dfPBC)[c(11,12)] <-c("Obs.16sC","Shan.16sC") #rename some columns

#Add alpha difference cols
dfPBC$ObsDiff<-dfPBC$Obs.16sW-dfPBC$Obs.16sC
dfPBC$ShanDiff<-dfPBC$Shan.16sW-dfPBC$Shan.16sC

#add columns for number of taxa gained and lost for each rat
#for each rat: count of taxa only in wild, count of taxa only in captive
output<-matrix(0, nrow =length(dfPBC$Woodrat.ID),  ncol =7 )
colnames(output)<-c("Woodrat.ID", "uniqueW", "lenW", "uniqueC", "lenC", "SharedWC", "TotalWC")

for(i in 1:length(dfPBC$Woodrat.ID)){
	rat<-as.character(dfPBC$Woodrat.ID[i]) #factors make it angry
	output[i,1]<-rat
	pw1<-subset_samples(psW, Woodrat.ID==rat)
	pw2<-prune_taxa(taxa_sums(pw1) > 0, pw1)
	taxW<-rownames(tax_table(pw2)) #all the taxa in the wild sample

	pc1<-subset_samples(psC, Woodrat.ID==rat)
	pc2<-prune_taxa(taxa_sums(pc1) > 0, pc1)
	taxC<-rownames(tax_table(pc2)) #all the taxa in the captive sample

	uniqueC<-length(setdiff(taxC,taxW))
	uniqueW<-length(setdiff(taxW,taxC))
	Shared<-length(intersect(taxW,taxC)) #what they have in common
	Total<-length(union(taxW,taxC)) #all taxa from a host wild an captive

	output[i,2]<-uniqueW
	output[i,3]<-length(taxW)
	output[i,4]<-uniqueC
	output[i,5]<-length(taxC)
	output[i,6]<-Shared
	output[i,7]<-Total
	}
dfWC<-data.frame(output, stringsAsFactors=FALSE)

dfWC$PropUC<-as.numeric(dfWC$uniqueC)/as.numeric(dfWC$lenC)
dfWC$PropUW<-as.numeric(dfWC$uniqueW)/as.numeric(dfWC$lenW)

#merge uniques onto dfPBC dataframe
dfWCU<-merge(dfPBC, dfWC, by = "Woodrat.ID")

#add columns with similarity (bc,j,wu,u) between wild-captive samples from each host
#a loop would have been neater, but this only repeats four times...

bray_WC<-phyloseq::distance(ps.rout, "bray")
matWC<-as.matrix(bray_WC)
matWC[lower.tri(matWC)] <- NA  #only need upper triangle, or will duplicate all the pairs
row.names(matWC) <- sample_data(ps.rout)$code
colnames(matWC) <- sample_data(ps.rout)$code
#reshape to a long format, need tidyr package
matWC_long<-gather(as.data.frame(matWC), key="code2", value= "distance")
ID1<-c(rep(rownames(matWC),dim(matWC)[1]))
matWC_long["code1"]<-ID1
matWC_long2<-na.omit(matWC_long)# get rid of the NA values from lower triangle
#add woodrat.ID for code1 and code2
codes<-data.frame(code=sample_data(ps.rout)$code, woodrat.ID=sample_data(ps.rout)$Woodrat.ID)
WC_long<-merge(matWC_long2, codes, by.x = "code1", by.y = "code")
WC_long2<-merge(WC_long, codes, by.x = "code2", by.y = "code")
WC_long2$samerat<-ifelse(WC_long2$woodrat.ID.x==WC_long2$woodrat.ID.y, "SR", "DR")
WC_long2$samesamp<-ifelse(WC_long2$code2==WC_long2$code1, "SS", "DS")
pairs_BC<-WC_long2[which(WC_long2$samerat=="SR" & WC_long2$samesamp=="DS"),]
pairs_BCs<-data.frame(Woodrat.ID=pairs_BC$woodrat.ID.x, bray=pairs_BC$distance)
#add it to the main dataframe
dfWCUB<-merge(dfWCU, pairs_BCs, by = "Woodrat.ID")

#repeat to add other dissimilarity metrics
#note that this overwrites some names from above
jac_WC<-phyloseq::distance(ps.rout, "jaccard")
matWC<-as.matrix(jac_WC)
matWC[lower.tri(matWC)] <- NA
row.names(matWC) <- sample_data(ps.rout)$code
colnames(matWC) <- sample_data(ps.rout)$code
matWC_long<-gather(as.data.frame(matWC), key="code2", value= "distance")
ID1<-c(rep(rownames(matWC),dim(matWC)[1]))
matWC_long["code1"]<-ID1
matWC_long2<-na.omit(matWC_long)
WC_long<-merge(matWC_long2, codes, by.x = "code1", by.y = "code")
WC_long2<-merge(WC_long, codes, by.x = "code2", by.y = "code")
WC_long2$samerat<-ifelse(WC_long2$woodrat.ID.x==WC_long2$woodrat.ID.y, "SR", "DR")
WC_long2$samesamp<-ifelse(WC_long2$code2==WC_long2$code1, "SS", "DS")
pairs_J<-WC_long2[which(WC_long2$samerat=="SR" & WC_long2$samesamp=="DS"),]
pairs_Js<-data.frame(Woodrat.ID=pairs_J$woodrat.ID.x, Jac=pairs_J$distance)
#add it to the main dataframe
dfWCUBJ<-merge(dfWCUB, pairs_Js, by = "Woodrat.ID")

Wuni_WC<-phyloseq::distance(ps.rout, "wunifrac")#this is slowish
matWC<-as.matrix(Wuni_WC)
matWC[lower.tri(matWC)] <- NA
row.names(matWC) <- sample_data(ps.rout)$code
colnames(matWC) <- sample_data(ps.rout)$code
matWC_long<-gather(as.data.frame(matWC), key="code2", value= "distance")
ID1<-c(rep(rownames(matWC),dim(matWC)[1]))
matWC_long["code1"]<-ID1
matWC_long2<-na.omit(matWC_long)
WC_long<-merge(matWC_long2, codes, by.x = "code1", by.y = "code")
WC_long2<-merge(WC_long, codes, by.x = "code2", by.y = "code")
WC_long2$samerat<-ifelse(WC_long2$woodrat.ID.x==WC_long2$woodrat.ID.y, "SR", "DR")
WC_long2$samesamp<-ifelse(WC_long2$code2==WC_long2$code1, "SS", "DS")
pairs_WU<-WC_long2[which(WC_long2$samerat=="SR" & WC_long2$samesamp=="DS"),]
pairs_WUs<-data.frame(Woodrat.ID=pairs_WU$woodrat.ID.x, Wunif=pairs_WU$distance)
#add it to the main dataframe
dfWCUBJW<-merge(dfWCUBJ, pairs_WUs, by = "Woodrat.ID")

Uni_WC<-phyloseq::distance(ps.rout, "unifrac")#this is slowish
matWC<-as.matrix(Uni_WC)
matWC[lower.tri(matWC)] <- NA
row.names(matWC) <- sample_data(ps.rout)$code
colnames(matWC) <- sample_data(ps.rout)$code
matWC_long<-gather(as.data.frame(matWC), key="code2", value= "distance")
ID1<-c(rep(rownames(matWC),dim(matWC)[1]))
matWC_long["code1"]<-ID1
matWC_long2<-na.omit(matWC_long)
WC_long<-merge(matWC_long2, codes, by.x = "code1", by.y = "code")
WC_long2<-merge(WC_long, codes, by.x = "code2", by.y = "code")
WC_long2$samerat<-ifelse(WC_long2$woodrat.ID.x==WC_long2$woodrat.ID.y, "SR", "DR")
WC_long2$samesamp<-ifelse(WC_long2$code2==WC_long2$code1, "SS", "DS")
pairs_U<-WC_long2[which(WC_long2$samerat=="SR" & WC_long2$samesamp=="DS"),]
pairs_Us<-data.frame(Woodrat.ID=pairs_U$woodrat.ID.x, Unif=pairs_U$distance)
#add it to the main dataframe
capdf<-merge(dfWCUBJW, pairs_Us, by = "Woodrat.ID") #give it a reasonable name

#add column for percent change observed diversity in captivity: (C-W)/W
capdf$Per.Change.Obs<-100*(capdf$Obs.16sC-capdf$Obs.16sW)/capdf$Obs.16sW

#adds on site code (eg L17) as new column, short codes used in manuscript
codes<-data.frame(code=site_codes$code2, Site_Sp=site_codes$Site_Sp)
capdf2<-merge(capdf, codes, by="Site_Sp")
saveRDS(capdf2,"capdf2_4Aug21.rds")

#################################################################################
#use data frame for analyses
#Are more diverse wild communities linked to more diverse captive communities?
ma<-lm(capdf$Obs.16sC~ capdf$Obs.16sW) #sig

#factors predicting change in diversity, composition
#Table S4, repeat model testing for each  metric
#run for responses: ObsDiff, ShanDiff, bray, Jac
#remove no-plant rats for nested model comparisons with diet
#confirmed that outliers do not alter results

cNA<-na.omit(capdf2)
m1<-lm(cNA$ObsDiff~ cNA$Species + cNA$Obs.plant + cNA$Obs.16sW)
m2<-lm(cNA$ObsDiff~ cNA$Species  + cNA$Obs.16sW)
anova(m1,m2) #no diet effects

#examine for all models
plot(fitted(m2), residuals(m2, type="pearson"), main= "Residuals vs Fitted",
                          xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
qqnorm(residuals(m2))

#return to dataset including animals with no diet data now that diet is dropped
m2<-lm(capdf2$ObsDiff~ capdf2$Species  + capdf2$Obs.16sW)
m3<-lm(capdf2$ObsDiff~  capdf2$Obs.16sW)
anova(m2,m3, test="F")  # compare models with F test, retain Species

#for plots of changing alpha diversity (Figure 5E)
#for plot, use % change in observed ASVs
rank<-data.frame(tapply(capdf2$Per.Change.Obs, capdf2$Species, mean))
names(rank)<-c("mean.obs")
rank$Species<-rownames(rank)
rank1<- rank[order(rank$mean.obs),]
level_order_dom <- as.vector(rank1$Species)

spec.col<-c("#F7C908","#1B0121", "#B72467", "#BFA8E2",
            "#4D04BF","#FF0000","#EF8630")

ggplot(data=capdf2, aes( x=factor(Species,level=level_order_dom),
                         y=Per.Change.Obs, fill=Species))+
	geom_boxplot(outlier.color="grey")+
	scale_fill_manual(values=spec.col)+
	geom_hline(yintercept=0, linetype = "dashed")+
	theme_classic()+
	theme(axis.text.x = element_text(angle = 90), legend.position="none")

#to show significant differences between groups on plot
m4<-lm(capdf2$Per.Change.Obs~  capdf2$Species)
pairs(emmeans(m4, ~Species))# posthoc

#quantify number of taxa gained and lost across sets of all wild, all captive
#checked- similar values whether using "rarefied" or "un-rarefied" ASVs
Wr<-prune_taxa(taxa_sums(psW) > 0, psW)
Cr<-prune_taxa(taxa_sums(psC) > 0, psC)
taxW<-rownames(tax_table(Wr))
taxC<-rownames(tax_table(Cr))

length(setdiff(taxW,taxC)) #unique to wilds
length(setdiff(taxC,taxW))  #unique to captives

length(setdiff(taxW,taxC))/length(taxW)  #prop unique all wilds
length(setdiff(taxC,taxW))/length(taxC)   #prop unique all caps

#compare to individuals
#proportion of ASVS found in only wild or captive for each individual animal
summary(capdf) #gives mean propUC, mean propUW
sd(capdf$PropUC)
sd(capdf$PropUW)
