
setwd("<directory_location>") #input working directory here
library(ggplot2)
library(phyloseq)
require(dplyr)
require(ggrepel)
library(reshape2)
library(betareg)

#files
#for phyloseq object:
taxatab<-read.csv("seqtab_trnL_nochim_5April20.Classified.csv") #from plant classifier
seqtab<-read.csv("seqtab_trnL_nochim_5April20.csv", row.names=1)
samdf <- read.csv("trnL_metadata_6April20.csv", header=TRUE)
#Isotope data
SIM<-read.csv("SIMMR_phylosymbiosis.csv", row.names=1) #for wild animals
SIM2<-read.csv("SIMMR_cactrial.csv") #for captive cactus feeding trial

###############################################################################
#create phyloseq object
#prep sample data. ***make sure row names match sequence table row names
#add a new column with the sample names used in the sequence table, set first colum to be row names
samdf1<-cbind(rownames(seqtab), samdf)
samdf2 <- data.frame(samdf1, row.names = 1)
taxa<-as.matrix(data.frame(taxatab,row.names=1))

#make phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(samdf2), tax_table(taxa))
saveRDS(ps, file = "trnL_ps_24April20.rds")
######################################################
#remove taxa not assigned to Kingdom Viridiplantae
psV<-subset_taxa(ps, Kingdom== "Viridiplantae")

#remove obvious bait contaminants(apple, peanut, soy, oats), Keep NAs
psb <- subset_taxa(psV,
	Genus!="Malus" &
	Genus!="Avena" &
	Genus!="Arachis" &
	 Genus!="Glycine"|
	is.na(Genus))

#remove samples that belong to other projects, hybrids, etc
psY<-subset_samples(psb, Diet.1=="Y")
#Subset to only wild samples
Wps<-subset_samples(psY, Sample_type=="wild")
#remove rare taxa (diet elements) from each sample, must be 1% of an animal's diet
##for each sample: set taxa < 1% to zero
Wp<- transform_sample_counts(Wps, function(x) replace(x, x<(.01*sum(x)),0) )
#remove anything no longer present in database
Wp1<-prune_taxa(taxa_sums(Wp) > 0, Wp)

#look at read counts:
names<-sample_names(Wp1)
sample_sum_df <- cbind(names,data.frame(sum = sample_sums(Wp1)))
sort <- sample_sum_df[order(sample_sum_df$sum),]
mean(sort$sum)
sd(sort$sum)

##################################
#psb is mostly wild data, but also includes samples from captive animals fed
#100% cactus diets, plot: cactus poorly amplifies in these animals
psD<-subset_samples(psb, Sample_type=="captive")
psC<-subset_samples(psD, Diet=="Cactus")
psC<-prune_taxa(taxa_sums(psC) > 1, psC)
plot_bar(psC, x="Woodrat.ID", fill="Family", title="3 week 100% cactus diet")

##############################################################################
#can we visualize lower read counts in cactus heavy samples?
RB<-subset_samples(Wp1, Run=="1") #restrict to single sequencing run for comparison

n2<-sample_data(RB)
n2$Names<- paste0( n2$Site_Sp, "-", n2$Woodrat.ID)
Reads <- data.frame(n2$Woodrat.ID, n2$Site_Sp, n2$Names, data.frame(sum = sample_sums(RB)))
high<-mean(Reads$sum)+sd(Reads$sum)
low<-mean(Reads$sum)-sd(Reads$sum)

ggplot(data=Reads, aes( x=n2.Site_Sp, y=sum))+
	geom_point(aes(color=n2.Site_Sp))+
	 geom_hline(yintercept=mean(Reads$sum), lty=2) +
	geom_hline(yintercept=c(low, high), lty=2, color="red") +
	theme_minimal()+
	theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "none")

#################################################################################
#Detecting missing cactus
#To estimate proportion of CAM/C4 plants in diets based on barcode data,
#first make csv of the wild samples with plant Family:genus, relative abundance
#merge plants at genus level
psg<-tax_glom(Wp1, taxrank="Genus", NArm=FALSE)
psgr  = transform_sample_counts(psg, function(x) x / sum(x) )#relative abundance
OTU1 = as(otu_table(psgr), "matrix")
#grab the rat id and site_sp cols from the sample data
dat<-data.frame(sample_data(psgr))
dat1<-dat[,c(4,9)]
OTU2<-merge(OTU1,dat1 , by=0)

OTU3 <- data.frame(OTU2[,-1], row.names = OTU2[,1]) #move rownames back to rownnames
Site_Sp<-OTU3$Site_Sp  #save for later
Woodrat.ID<-OTU3$Woodrat.ID #save for later

OTU3[,c(114:115)]<-NULL  #get rid of non-taxa
list1<-as.vector(colnames(OTU3))

Taxa1 = as(tax_table(psg), "matrix")
taxaW<-Taxa1[list1,]
taxaM<-data.frame(taxaW) #needs to be dataframe for paste to work
tmp <- paste(taxaM$Family, taxaM$Genus, sep=':')  #create "Family:genus" name for each taxa
names(OTU3) <- c( tmp) #replace sequences with Family:genus names
OTU4<-cbind(Site_Sp, Woodrat.ID, OTU3)#put sample info back: eg, Site_Sp.
write.csv(OTU4, "wild_rat_diets_overview24Apr.csv") #save csv
#for each (plant) family:genus determine whether it is CAM/C3/C4/Unknown.
#create new version of wild_rat_diets_overview24Apr.csv" just replacing OTU names
# with photosynthetic pathway (done in excel) this file is "wild_rat_diets_113Tax.csv"
#use this to calculate barcode-expected max and min amount of C4/CAM in each diet
C3<-read.csv("wild_rat_diets_113Tax.csv", check.names = FALSE)
C4<-data.frame(C3, row.names=1, check.names=FALSE)
Site_Sp<-C4$Site_Sp
Woodrat.ID<-C4$Woodrat.ID
C4[,c(1:2)]<-NULL  #get rid of non numeric reference info cols
C5<-t(rowsum(t(C4), group = colnames(C4)))  #transpose and use rowsums to merge the columns
C6<-data.frame(Site_Sp, Woodrat.ID,C5)

#make new data frame using C6: C3min, c3max, camc4min, camc4max
#excluding cactus from CAM values
C3min<-C6$C3
C3max<-C6$C3 + C6$U  #if all unknowns are C3
CAM4min<-C6$C4+C6$CAM  #if unknowns are all cam and C4
CAM4max<-C6$C4+C6$CAM+C6$U  #all cam and C4, plus assume unknowns are CAM/C4

C8<-data.frame(Site_Sp, C3min, C3max, CAM4min, CAM4max)
C9<-aggregate(C8, by=list(Site_Sp),FUN=mean, na.rm=TRUE)
C10 <- C9[order(C9$CAM4min),]
names1<- paste(C10$Group.1)

ggplot(data=C10, aes(x=Group.1, y=CAM4min, ymin=CAM4min, ymax=CAM4max)) +
        geom_pointrange() +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Population") +
		ylab("min and max non-cactus C4+CAM") +
		scale_x_discrete(limits=names1)+  #orders it by the names in the ordered column
        theme_bw()  # use a white background

#isotope mixing models give CAM/C4 plant proportions
#C6 has predicted amounts of CAM/C4 plants based on TRNL
#use both to determine which rats are missing cactus in trnl data
SIM<-read.csv("SIMMR_phylosymbiosis.csv", row.names=1)
C11<-data.frame(C6,row.names=2)

C12<-merge(SIM,C11, by=0) #default only keeps rows that were in both datasets

# as for pop level comparisons, use CAM, C4 without cactus sequences
C12$estCAM4min<-C12$C4+C12$CAM   #+C12$CAC
C12$estCAM4max<-C12$C4+C12$CAM+C12$U  #C12$CAC
CCC<-C12[,c(1,10, 16,17,7)]  #keep just the sumed DNAmin, DNAmax, and Iso values
colnames(CCC)<-c("ID", "Site_Sp", "DNA_CAM4_min","DNA_CAM4_max", "Iso_CAM4") #rename

glCAM4<-C12$global.C4.CAM.SIMMR.proportion
C13<-data.frame(C12$Site_Sp, C12$estCAM4min, C12$estCAM4max, glCAM4)
C14 <- melt(C13,
             variable.name = "CAM.C4",
             id.vars = "C12.Site_Sp")

###############################################################################
#use cactus feeding trial isotope data to model
#relationship between consumed cactus and fecal isotopes
SIM2<-read.csv("SIMMR_cactrial.csv")
SIM3<-SIM2[,c(1,2,6,8)] #keep needed columns
#reshape for plotting
SIM4 <- melt(SIM3,
             variable.name = "model",
             id.vars = c("ID", "cactus"))

#models fit with global C4/CAM estimates,
#explore link functions uses library(betareg)
logitCa <- betareg(value ~ cactus, data = SIM4,subset = model == "global.C4CAM", link="cauchit")# better shape than logit
logitCL <- betareg(value ~ cactus, data = SIM4,subset = model == "global.C4CAM", link="cloglog")# widest quantiles of links
logitL <- betareg(value ~ cactus, data = SIM4,subset = model == "global.C4CAM", link="log")# highest log liklihood

#quick look at two models, very similar
ggplot(data=SIM4, aes(x=cactus, y=value, color=model))+
	geom_point()+
	 geom_line(aes(y = predict(logitG, SIM4), colour = "global.C4CAM")) +
		geom_line(aes(y = predict(logitC, SIM4), colour = "cac.model")) +
	 geom_smooth(aes(fill=model),method = "lm", se=FALSE, linetype=4, size=.3)+
	 ylab("Isotopes: estimated C4/CAM") + xlab("Proportion cactus in diet")+
	theme_bw()

#prediction line with new data, add quantiles to prediction beta distribution w cloglog
fk<-data.frame(seq(0,1,.01))
colnames(fk)<-"cactus"
iso = predict(logitCL, fk, type="response") #predicted response
qs<-predict(logitCL,fk, type = "quantile", at = c(0.025,  0.975)) #quantiles
fk2<-cbind(fk,iso, qs)

#for a different regression model with log link
iso2 = predict(logitL, fk, type="response") #predicted response
qs2<-predict(logitL,fk, type = "quantile", at = c(0.025,  0.975)) #quantiles
fk3<-cbind(fk,iso2, qs2)

#for a different regression model with cauchit link
iso3 = predict(logitCa, fk, type="response") #predicted response
qs3<-predict(logitCa,fk, type = "quantile", at = c(0.025,  0.975)) #quantiles
fk4<-cbind(fk,iso3, qs3)

#plot with  2.5% and 97.5% quantiles of the predicted beta distribution
#using global mixing model
plot(iso~cactus, data=fk2, type="l", ylim=c(0, 1), main="cloglog (pseudo R2= 0.93)")
lines(q_0.025~cactus, data=fk2,  col="green")
lines(q_0.975~cactus, data=fk2,  col="green")
points(cac.model~cactus, data=SIM3, pch=16)

#How do real data fit on to this?
#isotope values=y, predicted barcode non-cac c4&cam=x
#support for missing cactus in diet if isotope results
#are outside of quantiles predicted from (max) barcode non-cactus levels
#show quantiles from 3 models
shapes<-c(21:25, 21:25, 21,22) #help differentiate points by using variety of shapes

ggplot()+
	geom_point(data=CCC, aes(x=DNA_CAM4_max, y=Iso_CAM4, color=Site_Sp,
						fill=Site_Sp,shape=Site_Sp), size=3)+
		scale_shape_manual(values=shapes)+
	#geom_line(data=fk2, aes(x=cactus,y=iso), color="black" )+
	geom_ribbon(data=fk2,aes(x=cactus, ymin=q_0.025, ymax=q_0.975), alpha=.2)+
	geom_ribbon(data=fk3,aes(x=cactus, ymin=q_0.025, ymax=q_0.975), alpha=.2)+
	geom_ribbon(data=fk4,aes(x=cactus, ymin=q_0.025, ymax=q_0.975), alpha=.2)+
		geom_text_repel(data=filter(CCC, Iso_CAM4>0.25), 	#add IDs to pts with iso vals>0.25
			aes(x=DNA_CAM4_max, y=Iso_CAM4,label=ID, color=Site_Sp))+
	labs(title = "cam/c4 estimates from wild fecal DNA v isotopes, w/ beta-reg quantiles diet trials",
			x="DNA Estimated CAM/C4 diet proportion", y="Isotope Estimated CAM/C4 diet proportion") +
	theme_bw()

#points outside of the prediction quantiles= strong evidence for cactus feeding
################################################################################
#adjusting for missing cactus, start with one pop, caspers bryanti
dna<-C11  #from above
dna$estCAM4min<-dna$C4+dna$CAM   #min non-cactus CAm/C4
dna$estCAM4max<-dna$C4+dna$CAM+dna$U  #max non catus C4/CAM
dna$CAM4avg<-(dna$estCAM4min+dna$estCAM4max)/2

Iso2<-SIM[,c(2,6)]
colnames(Iso2)[2] <-"SIMMR"  #rename column from "global.C4.CAM.SIMMR.proportion"

#using psg the ps object with plants merged to genus for all wild samples
#add in a new taxa for Family: Cactaceae, genus: cactus_adj
cac<-subset_taxa(psg, Family=="Cactaceae") #get taxonomy info for cactus
cac2<-as.data.frame(as(tax_table(cac), "matrix"))
cac3<-cac2[1,] #borrow opuntia taxonomy
cac3$Genus<-"cactus_iso"  # replace genus with cactus_iso name
row.names(cac3)

#add this new cactus as a row on full plant taxonomy
taxa2<-as.data.frame(as(tax_table(psg), "matrix"))
#bind it as last row on psg taxa table, note that DNA sequence has a 1 at the end
taxa3<-rbind(taxa2, cac3)
#taxa name "CTCCTTTTTTCAAAAAAAAAGAAAAAAATAAGGGTTCAGAAAGCAAGAATAAAAAAAAAAAG1" to use in the OTU table

#set up to modify cactus reads in the otu table
#samples are rows, plant taxa are columns.
#Add new cactus column to end, with above name
OTU2<-as.data.frame(as(otu_table(psg), "matrix"))

#export csv to determine which animals to adjust, based isotope v DNA mismatch
OTU3<-OTU2[,2:3] #just keep the 2 first columns
#grab the rat id and site_sp cols from the sample data
dat<-data.frame(sample_data(psg))
dat1<-dat[,c(4,9)]
OTU4<-merge(OTU3,dat1 , by=0, sort=FALSE)
OTU4[,2:3]<-NULL

write.csv(OTU4, "cac_adj.csv")
#go through by hand and designate cactus adjustment as Y/N
#N = no extra cactus needed, Y= some amount of extra cactus needed
#points outside of the prediction quantiles= strong evidence for cactus feeding
cacYN<-read.csv("cac_adj_YN.csv")

dna[,2:8]<-NULL  #keep the CAM4 average column
cac2<-merge(cacYN, dna, by.x=3, by.y=0, all.x=TRUE, sort=FALSE)
cac2[,6:7]<-NULL
cac3<-merge(cac2, Iso2, by.x=1, by.y=0, all.x=TRUE, sort=FALSE)

cac3$Cadj<-NA  #create column to fill in cactus values
cac3$Cadj<-ifelse(cac3$cac.adj== "N", 0, cac3$SIMMR - cac3$CAM4avg)  #fills in everything with complete data
cac3$missing<-ifelse(is.na(cac3$Cadj)==TRUE, "Y", "N")  #missing values

#pops where some individuals are missing isotope data;
#from lytle lep, casp_bry, CV-alb, [rio_alb]
#calculate population average isotopes from all isotopes run from that pop
#caspers bryanti
IsoCB<- SIM[ which(SIM$population=="Casper's" & SIM$species == "N. bryanti"), ]
avCB<-mean(IsoCB$global.C4.CAM.SIMMR.proportion, na.rm=TRUE)

#castle valley average, also use for neighboring Rio_alb
IsoCV<- SIM[ which(SIM$population=="Castlevalley" & SIM$species == "N. albigula"), ]
avCV<-mean(IsoCV$global.C4.CAM.SIMMR.proportion, na.rm=TRUE) #

IsoLR<- SIM[ which(SIM$population=="Lytleranch"), ] #lytle ranch
avLR<-mean(IsoLR$global.C4.CAM.SIMMR.proportion, na.rm=TRUE) #

#dumb way to do this...
mic<-cac3[which(cac3$missing=="Y"),]
mic$popav<-c(avLR,avCV,avCV, avCB, rep(avCV,7))
mic$Cadj<- mic$popav-mic$CAM4avg
mic[,11]<-NULL

cac4<-subset(cac3, !is.na(cac3$Cadj)) #remove the rows with NAs
cac5<-rbind(cac4, mic)
cac6<-cac5[,c(3,9)]
rownames(cac6) <- cac6[,1]

#Cadj gives the corrected total amount of cactus,
#but some animals already have some cactus reads, account for this
#for each animal need: total reads,  existing cactus reads
Trds<-data.frame(sample_sums(psg))
C<- subset_taxa(psg, Family=="Cactaceae")
rds<-cbind(Trds, sample_sums(C))
names(rds)<-c("tot.reads", "cac.reads")

cac7<-merge(rds, cac6, by=0, sort=FALSE)
cac7[,4]<-NULL
#calculate how many extra cactus reads to add on accounting for existing reads
cac7$cac.extra<-round((cac7$cac.reads-cac7$tot.reads*cac7$Cadj)/(cac7$Cadj-1))

#add this column onto OTU2, matching the new column added to the taxa table (taxa3)
OTU2$CTCCTTTTTTCAAAAAAAAAGAAAAAAATAAGGGTTCAGAAAGCAAGAATAAAAAAAAAAAG1<-cac7$cac.extra
taxa4<-as.matrix(taxa3)#taxa table must be a matrix

#put ps object back together with new taxa table and otu_table#################
psCacE <- phyloseq(otu_table(OTU2, taxa_are_rows=FALSE), sample_data(psg), tax_table(taxa4))
##for each sample: set taxa < 0 to zero,
#fixes theoretical risk that adjusted cactus is a negative amount
psCacE1<- transform_sample_counts(psCacE , function(x) replace(x, x<0,0) )
#remove two contaminated samples, previosly removed from 16s data
Samples_toRemoveP <- c("S83_1055w", "S84_1056w")
psCacE2<-subset_samples(psCacE1, !(code %in% Samples_toRemoveP))
#for each sample: set taxa < 1% to zero, done to original data,
#run to make sure added cactus pops are treated equivalently
psCacE3<- transform_sample_counts(psCacE2, function(x) replace(x, x<(.01*sum(x)),0) )
#remove anything no longer present in database
psCacE4<-prune_taxa(taxa_sums(psCacE3) > 0, psCacE3)
psCacE5<-rarefy_even_depth(psCacE4, 3600) #removes S10-955w, S58-1025w

saveRDS(psCacE5, "diet1percent_rar_29May20.rds")
