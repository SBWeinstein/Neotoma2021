#Test if low abundance taxa from wild are more likely to be lost in captivity

setwd()  # set file location

library(phyloseq)
library(tidyr)
library(glmmTMB)

#files
ps.rout<-readRDS( "ps.rout_phylo_18Oct2020.rds")

psC<-subset_samples(ps.rout, Sample_type=="captive")#subset to captive
psW<-subset_samples(ps.rout, Sample_type=="wild")#subset to wild

#only use animals w/ paired wild and captive 16s data
WC<-intersect(sample_data(psC)$Woodrat.ID, sample_data(psW)$Woodrat.ID) #125 rats
psWC<-subset_samples(ps.rout, (Woodrat.ID %in% WC))  #250 samples
psW125<-subset_samples(psWC, Sample_type=="wild")

#convert counts to relative abundance
psWr  = transform_sample_counts(psW125, function(x) x / sum(x) )

dfW<-data.frame(otu_table(psWr))
row.names(dfW)<-sample_data(psWr)$Woodrat.ID
dfW2<-data.frame(Woodrat.ID=sample_data(psWr)$Woodrat.ID , dfW)

#rename asvs with short codes A1......
asvs<-paste0("A",seq_len(dim(otu_table(psWr))[2]))
colnames(dfW2)<-c("Woodrat.ID", asvs)

#three columns: host ID, ASV, Relative abundance in wild sample
gW<-gather(dfW2, "ASV", "RA", -Woodrat.ID)

psC125<-subset_samples(psWC, Sample_type=="captive")

#for captive dataframe, transform to presence absence
psCp  = transform_sample_counts(psC125, function(x) ifelse(x>0,1,0) )

dfC<-data.frame(otu_table(psCp))
colnames(dfC)<-asvs

dfC2<-data.frame(Woodrat.ID=sample_data(psCp)$Woodrat.ID , dfC)
gC<-gather(dfC2, "ASV", "PA", -Woodrat.ID)
gCW<-merge(gW, gC, by.x=c("Woodrat.ID", "ASV"), by.y=c("Woodrat.ID", "ASV"))

#only interested in what wild asvs did.remove rows where the wild RA was 0
gCW2<-subset(gCW, gCW$RA!=0)

plot(gCW2$PA~log10(gCW2$RA))
hist(log10(gCW2$RA))

#model
dat.glmmTMBL <- glmmTMB(PA ~ log10(RA)+ (1|Woodrat.ID) ,
                        data = gCW2, family = "binomial")
