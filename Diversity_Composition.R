#based partially on:
#https://f1000research.com/articles/5-1492
#tutorials on https://joey711.github.io/phyloseq/index.html

library(ggplot2)
library(phyloseq)
library(vegan)
library(glmmTMB)
library(bbmle) ## for AICtab to compare models
library(ggeffects)# for plotting mixed effects model predictions #requires R version >3.4.3

setwd("<directory_location>") #input path here

#files
ps.rout<-readRDS("ps.rout_phylo_18Oct2020.rds") #rarefied, wild and captive rat 16s sequences
psCac<-readRDS("diet1percent_rar_29May20.rds") #diet

#subset to wild microbiome samples
psW<-subset_samples(ps.rout, Sample_type=="wild")
#create data frame with 16s diversity info
df <- data.frame(estimate_richness(psW, measures = c("Observed", "Shannon")),
	              Woodrat.ID=sample_data(psW)$Woodrat.ID,
	              Species=sample_data(psW)$Species,
	              site=sample_data(psW)$site,
	              Site_Sp=sample_data(psW)$Site_Sp,
	              conc=sample_data(psW)$ng.ul,
	              run=sample_data(psW)$Run)

#add columns for diet observed taxa, shannon, if specialist/generalist
psF<-tax_glom(psCac, "Family") #taxa glom plant data at family level
df_plant<-data.frame(estimate_richness(psF, measures=c("Observed", "Shannon")),
				Woodrat.ID=sample_data(psF)$Woodrat.ID,
				Sample_type=sample_data(psF)$Sample_type,
				Site_Sp=sample_data(psF)$Site_Sp)

#add a column for specialist designation,
#where single plant family >60% of diet (See Shipley 2009)
#merge each pop to get dominant plant, reads already equal between samples
merg = merge_samples(psF, "Site_Sp")
#convert counts to relative abundance
psFr = transform_sample_counts(merg, function(x) x/sum(x))

OT<-data.frame(otu_table(psFr))
#get the most abundant plant from each sample
plants<-colnames(OT)[max.col(OT,ties.method="first")]
#get max plant column index value (pick first if ties)
maxP<-max.col(OT, "first")
#get value from each of those columns
value <- OT[cbind(1:nrow(OT), maxP)]
res <- data.frame(Site_Sp=row.names(OT), plants, value)
#specialist defined as 60% of diet
res$SG60<-ifelse(res$value > 0.6, "S", "G")
#specialist defined as 50% of diet
res$SG50<-ifelse(res$value > 0.5, "S", "G")
#add spec/gen colum to plant dataframe
df_plant2<- merge(df_plant, res, by="Site_Sp")
#remove unneeded columns
df_plant2[c(1,5:7)]<-NULL

#merge plant diversity and 16s diversity into 1 dataframe,
#note that missing plant data for 8 16s samples
df16P<-merge(df, df_plant2, by= "Woodrat.ID", all.x=TRUE)
df16P<-merge(df, df_plant2, by= "Woodrat.ID")
no.plants<-df16P[which(is.na(df16P$SG50)==TRUE),]
no.plant.rats<-no.plants$Woodrat.ID
#make sequencing run a factor
df16P$run<-as.factor(df16P$run)

############################
#Do diet richness and Shannon diversity correlate? (using plant families)
cor.test(df16P$Observed.y, df16P$Shannon.y)

#Does diet richness vary among populations
mod1<-glm(df16P$Observed.y~df16P$Site_Sp + df16P$run,family=poisson, data=df16P)
mod2<-glm(df16P$Observed.y~df16P$Site_Sp,family=poisson, data=df16P)
anova(mod1, mod2, test="Chisq")# including run does not improve model
anova(mod2, test="Chisq") #populations sig vary

#Does diet composition vary among populations
#Permanova, adonis function, testing for different centroids
#repeat for bray and jaccard distance
psF_bray <- phyloseq::distance(psF, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psF))
#Adonis test
adonis(psF_bray ~ Site_Sp, data = sampledf)

############################################
#does diet diversity predict microbiome diversity?
#using DNA concentration as a proxy for DNA quality
fit<-glmmTMB(Observed.x~ Observed.y + conc +
              (1|Species)+ (1|Site_Sp) + (1|run),
               	data=df16P,
	             family="nbinom2")

#examine model fit, showing best model
plot(fitted(fit), residuals(fit, type="pearson"),
     main= "Residuals vs Fitted", xlab = "Fitted Values", ylab = "Residuals")
     abline(h = 0, lty = 2)
qqnorm(residuals(fit))
#compare AIC for models with different error distributions
AICtab(fit1, fit4)
#compare walds chisq to compare nested models with different variables
Anova(fit7, fit8)
#model outputs
summary(fit)

#Figure2: diet v microbiome diversity
#For regression line
pr <- ggpredict(fit, "Observed.y", type="random")
predicted<-data.frame(Observed.x=pr$x, Pred.y= pr$predicted)

spec.col<-c("#F7C908","#1B0121", "#B72467", "#BFA8E2", "#4D04BF","#FF0000","#EF8630")

ggplot()+
	scale_color_manual(values = spec.col)+
	scale_shape_manual(values = c(16,17,17,18,19,18,15))+
	geom_point(data=df16P, aes(y=Observed.x, x= Observed.y, color=Species, shape=Species),
	                            size=3, alpha=.8)+
	geom_line(data=predicted,aes(y=Pred.y, x= Observed.x), size=1)+
	xlab("Diet diversity") + ylab("Microbiome richness")+
	xlim(1,13)	+
	theme_classic()	+
	theme(legend.position = "none")

#####################
##does microbiome dispersion increases with diet diversity?
#make a data frame from the sample_data
#sampledf <- data.frame(sample_data(psW))
psW_bray<-phyloseq::distance(psW, "bray")
psW_wu<-phyloseq::distance(psW, "wunifrac")
psW_jac<-phyloseq::distance(psW, "jaccard")
psW_u<-phyloseq::distance(psW, "unifrac")

#calculate dispersion
bm1<-betadisper(psW_bray, sampledf$Site_Sp, bias.adjust = TRUE)
bm2<-betadisper(psW_jac, sampledf$Site_Sp, bias.adjust = TRUE)
bm3<-betadisper(psW_wu, sampledf$Site_Sp, bias.adjust = TRUE)
bm4<-betadisper(psW_u, sampledf$Site_Sp, bias.adjust = TRUE)

#get average distance to group, repeat for each distance metric
mod1<-data.frame(with(bm4, tapply(distances, sampledf$Site_Sp, "mean")))
obsD<-data.frame(with(df16P, tapply( Observed.y, Site_Sp, "mean")))
colnames(obsD)<-"meanDO"
merg<-merge(mod1, obsD, by=0)
#add species to datframe.
spec<-data.frame(df16P$Site_Sp, df16P$Species)
spec2<-unique(spec)
colnames(spec2)<-c("Site_Sp", "Species")
merg1<-merge(merg, spec2, by=1)

fit3 <- glmmTMB(disBray~ meanDO + (1|Species),
	          data=merg1,
	          family=beta_family())
#examine model fits, compare models as above
