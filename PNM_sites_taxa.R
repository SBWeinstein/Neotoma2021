#test whether increased host species diversity at a site reduces
#fit of Prokaryote Neutral Model
#code for Prokaryote Neutral Model from the supplemental material of 
#Burns, A., Stephens, W., Stagaman, K. et al. Contribution of neutral processes 
#to the assembly of gut microbial communities in the zebrafish over host development. 
#ISME J 10, 655–664 (2016). https://doi.org/10.1038/ismej.2015.142


setwd()

library(ggplot2)
library(phyloseq)
library(plyr) #count function
library(minpack.lm) #for sncm.fit function 
library(Hmisc) #for sncm.fit function
library(stats4)#for sncm.fit function

#files
ps.rout<-readRDS("ps.rout_phylo_18Oct2020.rds")

###sncm.fit function from Burns et al. 
#############################################################################
sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
	options(warn=-1)

	#Calculate the number of individuals per community
	N <- mean(apply(spp, 1, sum))
	
	#Calculate the average relative abundance of each taxa across communities
	if(is.null(pool)){
		p.m <- apply(spp, 2, mean)
		p.m <- p.m[p.m != 0]
		p <- p.m/N
	} else {
		p.m <- apply(pool, 2, mean)
		p.m <- p.m[p.m != 0]
		p <- p.m/N
	}

	#Calculate the occurrence frequency of each taxa across communities
	spp.bi <- 1*(spp>0)
	freq <- apply(spp.bi, 2, mean)
	freq <- freq[freq != 0]

	#Combine
	C <- merge(p, freq, by=0)
	C <- C[order(C[,2]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
	p <- C.0[,2]
	freq <- C.0[,3]
	names(p) <- C.0[,1]
	names(freq) <- C.0[,1]

	#Calculate the limit of detection
	d = 1/N

	##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
	m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
	m.ci <- confint(m.fit, 'm', level=0.95)
	
	##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
	sncm.LL <- function(m, sigma){
		R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
		R = dnorm(R, 0, sigma)
		-sum(log(R))
	}
	m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
	
	##Calculate Akaike's Information Criterion (AIC)
	aic.fit <- AIC(m.mle, k=2)
	bic.fit <- BIC(m.mle)

	##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
	freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
	Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
	
	pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for binomial model
	bino.LL <- function(mu, sigma){
		R = freq - pbinom(d, N, p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.bino <- AIC(bino.mle, k=2)
	bic.bino <- BIC(bino.mle)
	
	##Goodness of fit for binomial model
	bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
	Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))

	bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for Poisson model
	pois.LL <- function(mu, sigma){
		R = freq - ppois(d, N*p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.pois <- AIC(pois.mle, k=2)
	bic.pois <- BIC(pois.mle)
	
	##Goodness of fit for Poisson model
	pois.pred <- ppois(d, N*p, lower.tail=FALSE)
	Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))

	pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

	##Results
	if(stats==TRUE){
		fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
		fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
		return(fitstats)
	} else {
		A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
		A <- as.data.frame(A)
		colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
		if(is.null(taxon)){
			B <- A[order(A[,1]),]
		} else {
			B <- merge(A, taxon, by=0, all=TRUE)
			row.names(B) <- B[,1]
			B <- B[,-1]
			B <- B[order(B[,1]),]
		}
		return(B)
	}
}

################################################################
#subset to wild only
psW<-subset_samples(ps.rout, Sample_type=="wild")
psW<-prune_taxa(taxa_sums(psW) > 0, psW) #2704 in all wilds

samdf<-sample_data(psW)

#get sites for analysis
samdf$site[samdf$site == "Pioneertown Reserve"] <- "Pioneertown"   #groups these two pops for this analysis only
sample_data(psW) <- samdf  #replace dataframe with altered site names
sites<-count(samdf, "site", ) #how many animals sampled at each site?  
sites8<-sites[sites$freq > 6,] #this gives 15 sites with  7 animals
s8_list<-sites8$site  #list of these sites

#loop thru all the sites to get table of neutral model R2 by site
output<-data.frame(0, nrow=length(s8_list), ncol=4)

for(i in 1:length(s8_list)){
	pop<-s8_list[i]
	output[i,1]<-as.character(pop)
	pspop<-subset_samples(psW, site==pop) #using sites
	pspop<-prune_taxa(taxa_sums(pspop) > 0, pspop) 
	spp<-data.frame(otu_table(pspop))
	model<-sncm.fit(spp=spp, pool=NULL, stats=TRUE, taxon=NULL)
	
	output[i,2]<-nsamples(pspop)
	output[i,3]<-length(unique(sample_data(pspop)$Species))
	output[i,4]<-model$Rsqr
	}
colnames(output)<-c("site", "nrats", "nspp" ,"R2")

summary(lm(output$R2~output$nspp)) #more species at a site, worse R2

#plot for figure S5B
ggplot(data=output, aes(x=nspp, y=R2))+
	geom_smooth(method='lm', formula= y~x)+	
	geom_point(aes(size = nrats), alpha=0.8)+
	scale_size(range = c(7/2, 16/2))

#compare captive v wild PNM fit
spp<-data.frame(otu_table(psW)) #all ASVs, switch to psC for captive R2
#model fit
test<-sncm.fit(spp=spp, pool=NULL, stats=TRUE, taxon=NULL)

#taxa info, for making plots, and pulling out neutral/non taxa
test3<-sncm.fit(spp=spp, pool=NULL, stats=FALSE, taxon=NULL)

#which taxa are above, below, within predictions
test3$low<- ifelse(test3$freq<test3$pred.lwr, 0, 1)
test3$high<-ifelse(test3$freq>test3$pred.upr, 1, 0)
test3$vals<-test3$high + test3$low

test3$vals[test3$vals == 0] <- "Below"   
test3$vals[test3$vals == 1] <- "Within"   
test3$vals[test3$vals == 2] <- "Above"  
test3$vals<-factor(test3$vals, levels =c("Above" ,  "Within" ,"Below"))

###pull out list of ASVs that are neutral and non-neutral, above and below
#these can be used to subset taxa for MRM models
all<-rownames(test3) #2704 taxa
neutral<-rownames(test3[test3$vals=="Within",]) 
non_neu<-setdiff(all,neutral) 

#pull out "selected for/against" taxa
Above<-rownames(test3[test3$vals=="Above",]) 
Below<-rownames(test3[test3$vals=="Below",]) 
