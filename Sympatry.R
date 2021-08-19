#does sympatry increases ASV sharing?
#compare pair-wise Jaccard dissimilarities for 
#heterospecifics at the same site versus equivalent pairings at different sites
#then compare average percent asvs shared between sympatric populations
#to average shared between equivalent allopatric pops

setwd("<>") #set location of files

library(tidyr)
library(phyloseq)
library(ggplot2) # for plot

#files
ps.rout<-readRDS("ps.rout_phylo_18Oct2020.rds")
codes<-read.csv("Site_Sp_codes.csv")

#subset to wild
psW<-subset_samples(ps.rout, Sample_type=="wild")

#similarity based on ASV presence
jac_W<-phyloseq::distance(psW, "jaccard")#All possible combos
matW<-as.matrix(jac_W)
matW[lower.tri(matW)] <- NA  #only need upper triangle, or will duplicate all the pairs

row.names(matW) <- sample_data(psW)$code
colnames(matW) <- sample_data(psW)$code

matW_long<-gather(as.data.frame(matW), key="code2", value= "distance")
ID1<-c(rep(rownames(matW),dim(matW)[1]))
matW_long["code1"]<-ID1
matW_long2<-na.omit(matW_long)# get rid of the NA values from lower triangle

#add woodrat.ID for code1 and code2
codes<-data.frame(code=sample_data(psW)$code, woodrat.ID=sample_data(psW)$Woodrat.ID)
W_long<-merge(matW_long2, codes, by.x = "code1", by.y = "code")
W_long2<-merge(W_long, codes, by.x = "code2", by.y = "code")

#add in columns with the sites and species for each sample
mdat<-data.frame(code=sample_data(psW)$code, Species=sample_data(psW)$Species, site=sample_data(psW)$site)

W3<-merge(W_long2, mdat, by.x = "code1", by.y = "code")#add the site and spp info for rat code1

W4<-merge(W3, mdat, by.x = "code2", by.y = "code")#add the site and spp info for rat code2

#clean up column names
colnames(W4)<-c("code2", "code1", "distance", "rat1", "rat2", "sp1", "site1", "sp2", "site2")

#set up columns for subsetting for analyses
W4$SameRat<-ifelse(W4$rat1==W4$rat2, "SR", "DR") #SR indicates its the same sample/Rat (BC==0)
W4$SameSite<-ifelse(W4$site1==W4$site2, "SS", "DS") #from the same site?
W4$SameSpp<-ifelse(W4$sp1==W4$sp2, "SS", "DS") #from the same species?

#replicates are limited to bryanti, lepida, and macrotis, restrict dataframe to those
blm<-c("N. lepida", "N. macrotis", "N. bryanti")

W5<-W4[which(W4$sp1 %in% blm & W4$sp2 %in% blm),]
W6<-W5[which(W5$SameRat == "DR"),] #remove the self comparisons

#caculate similarities for species pairs that are sympatric and not 
#B and L
#sympatric B and L (site is White Water)
WW<-W6[which(W6$site1== "White Water" & W6$SameSite=="SS" & W6$SameSpp=="DS"),]
WWBL<-WW$distance  #average distance between symptratic L and B at WW

#average distance between B and L, allopatric 
bl<-c("N. lepida",  "N. bryanti")
BL<-W6[which(W6$sp1 %in% bl & W6$sp2 %in% bl & W6$SameSite=="DS" & W6$SameSpp=="DS"),]
WBL<-BL$distance

t.test(WBL,WWBL)
#Jaccard: t = 10.806, df = 38.396, p-value = 3.324e-13

#store to make plot S4
dfplot<-data.frame(distance=WWBL, pair=rep("BL", length(WWBL)), Sym=rep("Sym", length(WWBL)))
ABL<-data.frame(distance=WBL, pair=rep("BL", length(WBL)), Sym=rep("Allo", length(WBL)))
dfplot1<-rbind(dfplot, ABL)

#repeat for caspers macrotis and bryanti
ca<-W6[which(W6$site1== "Casper's" & W6$SameSite=="SS" & W6$SameSpp=="DS"),]
caMB<-ca$distance  #average distance between symptratic M and B at WW

#average distance between B and M, allopatric 
mb<-c("N. macrotis",  "N. bryanti")
MB<-W6[which(W6$sp1 %in% mb & W6$sp2 %in% mb & W6$SameSite=="DS" & W6$SameSpp=="DS"),]
aloMB<-MB$distance

t.test(caMB,aloMB) #Bray: t = -2.6858, df = 31.2, p-value = 0.01149
#Jaccard: t = -2.638, df = 31.16, p-value = 0.0129

#setting up plot, continued
dfp1<-data.frame(distance=caMB, pair=rep("MB", length(caMB)), Sym=rep("Sym", length(caMB)))
dfp2<-data.frame(distance=aloMB, pair=rep("MB", length(aloMB)), Sym=rep("Allo", length(aloMB)))
dfplot2<-rbind(dfplot1, dfp1, dfp2)

#repeat for big morongo macrotis and lepida
Mor<-W6[which(W6$site1== "Big Morongo" & W6$SameSite=="SS" & W6$SameSpp=="DS"),]
MorML<-ca$distance  #average distance between symptratic L and B at WW

#average distance between L and M, allopatric 
MacL<-c("N. macrotis",  "N. lepida")
ML<-W6[which(W6$sp1 %in% MacL & W6$sp2 %in% MacL & W6$SameSite=="DS" & W6$SameSpp=="DS"),]
aloML<-ML$distance

t.test(MorML,aloML) #Bray: t = -2.632, df = 31.031, p-value = 0.01311
#jaccard: t = -2.6048, df = 31.021, p-value = 0.01399

dfp3<-data.frame(distance=MorML, pair=rep("ML", length(MorML)), Sym=rep("Sym", length(MorML)))
dfp4<-data.frame(distance=aloML, pair=rep("ML", length(aloML)), Sym=rep("Allo", length(aloML)))
dfplot3<-rbind(dfplot2, dfp3, dfp4)

#make plot
theme_set(theme_bw())
labels <- c(BL = "bryanti-lepida", MB = "macrotis-bryanti", ML="macrotis-lepida")

ggplot(dfplot3, aes(x=Sym, y=distance)) + 
	geom_jitter(alpha=c(.08)) +
	stat_summary(fun.data=mean_sdl, aes(color=Sym))+
	facet_grid(.~pair, labeller=labeller(pair = labels))+
	ylab("Jaccard Distance")+
	theme( legend.position = "none", axis.title.x = element_blank())

############################################################################
#calculate average percent asvs shared between sympatric populations at each site
#then calculate average shared between equivalent allopatric pops

#subset to only B,M,L species
LMB<-c("N. lepida", "N. macrotis", "N. bryanti")
psLMB<-subset_samples(psW, Species %in% LMB)
psLMB<-prune_taxa(taxa_sums(psLMB) > 0, psLMB)

#merge samples by population
LMB_pops = merge_samples(psLMB, "Site_Sp") #rownames are Site_Sp, other sample_data garbage
sample_data(LMB_pops)$Site_Sp<-rownames(sample_data(LMB_pops)) #put Site_Sp back in
pops<-c(rownames(sample_data(LMB_pops)))

#get the taxa in each  population
pops[1]
outlist<-list()
for(i in 1:length(pops)){
	pop1<-pops[i]
	pw1<-subset_samples(LMB_pops, Site_Sp==pop1)
		pw2<-prune_taxa(taxa_sums(pw1) > 0, pw1)
		taxW<-rownames(tax_table(pw2)) #all the taxa in pop1
	outlist[[i]]<-taxW
	}

########get overlap between each pop
output<-matrix(0, nrow =0,  ncol =7 )
colnames(output)<-c("pop1", "pop2", "ASV1", "ASV2", "shared", "to", "pshare")

for(j in 1:length(pops)){

	taxA<-outlist[[j]]
	output2<-matrix(0, nrow =19,  ncol =7 )
	colnames(output2)<-c("pop1", "pop2", "ASV1", "ASV2", "shared", "to", "pshare")
	
	for(i in 1:length(pops)){
		taxB<-outlist[[i]]
		Asv1<-length(taxA)
		Asv2<-length(taxB)
		shared<-length(intersect(taxA,taxB))#how many are shared
		tot<-length(union(taxA,taxB))#how many there are at the site
		pshare<-shared/tot
	
		output2[i,1]<-pops[j]  #pop1 (new)
		output2[i,2]<-pops[i]	#pop2
		output2[i,3]<-Asv1	#number of asvs in 1
		output2[i,4]<-Asv2
		output2[i,5]<-shared
		output2[i,6]<-tot
		output2[i,7]<-pshare
		}
	output<-rbind(output,output2)
	}

dfpops<-data.frame(output, stringsAsFactors=FALSE)
#dfpops$p21<-paste(dfpops$pop2, dfpops$pop1)
#hacky, but we expect 171 combinations of 2 picked from 19
#length(unique(dfpops$shared))

#need to get rid of the replicates#make it into a matrix and get rid of lower triangle?
#library(tidyr)
df2<-dfpops[,c(1,2,7)]

#add on species1, species2, site1, site2
#codes<-read.csv("Site_Sp_codes.csv")
cd1<-codes[,c(1,4,10)]
df9<-merge(df2, cd1, by.x="pop1", by.y="Site_Sp")
df10<-merge(df9, cd1, by.x="pop2", by.y="Site_Sp")

df10$SSp<-ifelse(df10$Species.x==df10$Species.y, "SR", "DR")#same species?
df10$SSt<-ifelse(df10$Location.x==df10$Location.y, "SS", "DS")#same site?

lep<-df10[which(df10$Species.x=="N. lepida"),]  #pop1 is lepida

#white water LB
lepbry<-lep[which(lep$Species.y=="N. bryanti"),] #pop2 is bryanti
lepbryDS<-lepbry[which(lepbry$SSt=="DS"),] #get the ones from diff sites
lepbrySS<-lepbry[which(lepbry$SSt=="SS"),] #get the ones from diff sites

#allopatric
mean(as.numeric(lepbryDS$pshare))#0.2598178
sd(as.numeric(lepbryDS$pshare))#0.06708577

#Sympatric
mean(as.numeric(lepbrySS$pshare))#0.3940193, a 51.7% increase

#BM lep v mac
lepmac<-lep[which(lep$Species.y=="N. macrotis"),] #pop2 is N. macrotis
lepmacDS<-lepmac[which(lepmac$SSt=="DS"),] #get the ones from diff sites

mean(as.numeric(lepmacDS$pshare))#0.1045157
sd(as.numeric(lepmacDS$pshare))#0.01495982

lepmacSS<-lepmac[which(lepmac$SSt=="SS"),] #get the ones from same site
mean(as.numeric(lepmacSS$pshare))#0.093558282208589 shared fewer than average

#casp bry v mac
bry<-df10[which(df10$Species.x=="N. bryanti"),]  #pop1 is N. bryanti

brymac<-bry[which(bry$Species.y=="N. macrotis"),] #pop2 is N. macrotis
brymacDS<-brymac[which(brymac$SSt=="DS"),] #get the ones from diff sites

mean(as.numeric(brymacDS$pshare))# 0.1242243
sd(as.numeric(brymacDS$pshare))#0.04066798

brymacSS<-brymac[which(brymac$SSt=="SS"),] #get the ones from same site
mean(as.numeric(brymacSS$pshare))# 0.280874316939891 more than average, increase of 133.3%



