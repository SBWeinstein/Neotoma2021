#phylogeny and geography contributions
#examined using (partial) mantel tests and mutual clustering approaches
setwd("<>") #set working directory
library("vegan")
library("geosphere")
library("ecodist")
library("tidyr")
library("ggplot2")
library("phyloseq")
library("phangorn")
library("phylogram")
library("TreeDist")
library("TreeTools")

#files
ps.rout<-readRDS("ps.rout_phylo_18Oct2020.rds") # microbiome
site<-as.data.frame(read.csv("Site_lat_long.csv")) #lat/long data
Neo_tree<-read.tree(file="Neotoma_phylosymtree.nwk") #Neotoma tree with branch lengths from Matocq 2007
Neo_tree1<-read.tree(file="Neo_bif_tree2.phy") #rat tree for clustering analyses
sites<-as.data.frame(read.csv("Site_Sp_codes.csv", header=TRUE)) #metadata for populations

#subset to wild microbiome samples
psW<-subset_samples(ps.rout, Sample_type=="wild")

#Does composition differ by sites?
sampledf <- data.frame(sample_data(psW))
psW_bray<-phyloseq::distance(psW, "bray")
psW_wu<-phyloseq::distance(psW, "wunifrac")
psW_jac<-phyloseq::distance(psW, "jaccard")
psW_u<-phyloseq::distance(psW, "unifrac")

# Adonis test - site first, controlling for species
adonis(psW_bray ~ site+Species, data = sampledf)  #sig
adonis(psW_jac ~ site+Species, data = sampledf)  #sig
adonis(psW_wu ~ site+Species, data = sampledf)  #sig
adonis(psW_u ~ site+Species, data = sampledf)  #sig

# Adonis test - species first, controlling for site
adonis(psW_bray ~ Species+site, data = sampledf)  #sig
adonis(psW_jac ~ Species+site, data = sampledf)  #sig
adonis(psW_wu ~ Species+site, data = sampledf)  #sig
adonis(psW_u ~ Species+site, data = sampledf)  #sig

#########################
#(partial) mantel tests at the individual, pop and species levels
#create distance matrices
#make microbiome distance matrix
psW_bray<-phyloseq::distance(psW, "bray")  #use for full ASV set
psW_jac<-phyloseq::distance(psW, "jaccard") #use for subset comparisons
dm<-as.matrix(psW_bray) #155*155  # options: psW_jac, bray
dm2<-as.dist(dm)

#for geography
D_site<-distm(cbind(site$longitude,site$latitude), fun=distGeo) #longitude first
row.names(D_site) <- site$site
colnames(D_site) <- site$site
D_site_km<-D_site/1000 #convert to km
#convert to long format
D_site_long<-gather(as.data.frame(D_site_km), key="Site2", value= "distance")
D_site_long["Site1"]<-c(rep(rownames(D_site_km), 19))
D_site_long <- D_site_long[c(3,1,2)]
#for goegraphy make blank new matrix with dims and names matching the microbiome matrix
big.site.matrix <- matrix(0, nrow = length(ps_data$site), ncol = length(ps_data$site))
row.names(big.site.matrix) <- ps_data$site
colnames(big.site.matrix) <- ps_data$site
#use the table to fill in values in geo matrix with dimensions matching all the samples
for (i in 1:nrow(D_site_long)) {
   big.site.matrix[rownames(big.site.matrix) == D_site_long$Site2[i], colnames(big.site.matrix) == D_site_long$Site1[i]] <- D_site_long$distance[i]
}
gm<-as.dist(big.site.matrix)

#make distance matrix for phylogeny
#replace tip labels to march species names in other files
#note that as of oct20, extra space after stephensi fixed
Neo_tree$tip.label<-c("N. albigula", "N. stephensi", "N. macrotis",
			"N. cinerea", "N. bryanti", "N. lepida", "N. devia" )
D_evo<-cophenetic(Neo_tree)
#convert into long format, species1, species2, distance
D_evo_long<-gather(as.data.frame(D_evo), key="Species2", value= "distance")
D_evo_long["Species1"]<-c(rep(rownames(D_evo), 7)) #make extra column for compared to species
D_evo_long <- D_evo_long[c(3,1,2)]  #reorder columns
#for phylogeny make blank new matrix with dims and names matching the microbiome matrix
evo.matrix <- matrix(0, nrow = length(ps_data$site), ncol = length(ps_data$site))
row.names(evo.matrix) <- ps_data$Species
colnames(evo.matrix) <- ps_data$Species

#use the table to fill in values in phylo matrix with dimensions matching all the samples
for (i in 1:nrow(D_evo_long)) {
   evo.matrix[rownames(evo.matrix) == D_evo_long$Species2[i],
	colnames(evo.matrix) == D_evo_long$Species1[i]] <- D_evo_long$distance[i]
	}
dist.evo = as.dist(evo.matrix)

#partial mantel test, individual microbiome and geography,conditioned on phylogeny
pm_MG = vegan::mantel.partial(dm2, gm, dist.evo, method = "pearson", permutations = 9999)
#partial mantel test, individual microbiome and phylo, conditioned on geo
pm_MP = vegan::mantel.partial(dm2, dist.evo, gm, method = "pearson", permutations = 9999)

#diet, mantel test and partial mantel
diet_mant  = vegan::mantel(dietmat,dist.evo, method = "pearson", permutations = 9999, na.rm = TRUE)
diet_pm = vegan::mantel.partial(dm2, dietmat, dist.evo, method = "pearson", permutations = 9999)

##############################################################################
#are these patterns driven by really strong population differences?
#merge animals from each pop together
ps_pop<-merge_samples(psW, "Site_Sp")
ppr  = transform_sample_counts(ps_pop, function(x) x / sum(x) )#convert to relative abundance
ppr_dist<-phyloseq::distance(ppr, "bray")  # or bray,jaccard
dmp<-as.matrix(ppr_dist)
dmp2<-as.dist(dmp)

#replace the metadata, lost when merged
#sites<-as.data.frame(read.csv("Site_Sp_codes.csv", header=TRUE))
p_dat<-as.data.frame(sample_data(ppr))
p_dat1<-merge(p_dat, sites, by.x=0, by.y=10)

rownames(p_dat1)<-p_dat1$Row.names
geo2 = data.frame(p_dat1$Longitude, p_dat1$Latitude)
d.geo2 = distm(geo2, fun = distGeo)
dist.geo2 = as.dist(d.geo2)

#use partial mantel test with rat phylogeny
#make third distance matrix for rat species
# start with D_evo_long, from above
# rat species are in p_dat1$Species.y
#make blank new matrix with dims and names matching the microbiome matrix
evo.matrix <- matrix(0, nrow = length(p_dat1$Species.y), ncol = length(p_dat1$Species.y))
row.names(evo.matrix) <- p_dat1$Species.y
colnames(evo.matrix) <- p_dat1$Species.y

#use the table to fill in values in matrix with dimensions matching all the samples
for (i in 1:nrow(D_evo_long)) {
   evo.matrix[rownames(evo.matrix) == D_evo_long$Species2[i],
	colnames(evo.matrix) == D_evo_long$Species1[i]] <- D_evo_long$distance[i]
	}
dist.evo = as.dist(evo.matrix)

pop_gE  = vegan::mantel.partial(dmp2, dist.geo2, dist.evo, method = "pearson", permutations = 9999)
pop_p  = vegan::mantel.partial(dmp2, dist.evo, dist.geo2, method = "pearson", permutations = 9999)

#just looking at phylogeny, at pop level, no geographic control
pop_phy  = vegan::mantel(dmp2, dist.evo, method = "pearson", permutations = 9999, na.rm = TRUE)

#test phylogeny effects when merged at species level.
#merging sites, diets so now using Mantel, not partial mantel test
#need order of D_evo and dmp to match, do by hand, probably a smarter way...
D_ev2<-D_evo[c(1,5,4,7,6,3,2),c(1,5,4,7,6,3,2)]
D_ev3<-as.dist(D_ev2)

ps_sp<-merge_samples(psW, "Species")
spr  = transform_sample_counts(ps_sp, function(x) x / sum(x) )#convert to relative abundance
spr_dist<-phyloseq::distance(spr, "unifrac")# repeat for bray, jaccard, wunifrac, unifrac
smp<-as.matrix(spr_dist)
smp2<-as.dist(smp)

spp_phy  = vegan::mantel(smp2, D_ev3, method = "pearson", permutations = 999, na.rm = FALSE)
###############################################################################
###############################################################################
####Alternative to Mantel tests: test tree congruence
#test for significance using permutation test

#run first for host species
UP <- hclust(dmp2, method="average")  #build UPGMA tree
micros <- as.dendrogram(UP)
mic_new<-write.dendrogram(micros, edges = TRUE)
mic_tree<-ape::read.tree(text=mic_new)
#make sure tip names match
mic_tree$tip.label<-c("N. macrotis", "N. stephensi","N. cinerea","N. bryanti",
                      "N. lepida","N. albigula","N. devia" )

CID <- ClusteringInfoDistance(mic_tree, Neo_tree, normalize = FALSE)
ExpectedVariation(mic_tree, Neo_tree, samples = 100000)
distance <- TreeDistance(mic_tree, Neo_tree)

#using TreeTools package to bootstrap p-values
perm<-100000
output<-matrix(0, nrow = perm,  ncol =1)
for(i in 1:perm){
	RT1<-RandomTree(tips=7, root = FALSE)
	RT2<-RandomTree(tips=7, root = FALSE)
	CID_RT <- TreeDistance(RT1,RT2)
	output[i,1]<-CID_RT
	}

mean(output[,1])
#compute p-value
#count the number of values (statistics) that are greater than or equal
#to the observed value, and divide by the number of values
boots<-output[,1]
length(boots[boots> distance])/perm

sb<-sort(boots)
sb[(.025*perm)] #would be 90% CI with 5% on either ends, use 2.5% for 95%
sb[(.975*perm)]

#run for host populations
#requires bifurcating trees.
#use tree with placeholders of 2 pops per species, for species with >1 pop
#Neo_tree1<-read.tree(file="Neo_bif_tree2.phy")
ps_pop<-merge_samples(psW, "Site_Sp")
ppr  = transform_sample_counts(ps_pop, function(x) x / sum(x) )

#for merged populations, replace merged metadata to get site info
p_dat<-as.data.frame(sample_data(ppr))
p_dat1<-merge(p_dat, sites, by.x=0, by.y=10)#for merged pops
p_dat2<-data.frame(row.names=p_dat1$Row.names, Species=p_dat1$Species.y,
	                  pop=p_dat1$Row.names, code= p_dat1$code2)
sample_data(ppr) <- p_dat2

permut<-100000
dist_output<-matrix(0, nrow = permut,  ncol =1)
for(i in 1:permut){
	#Generalized RF methods need bifurcating trees, subsample pops to produce
	#tree with no more than two pops per species, matching Neo_bif_tree2
	#randomly pick 2 pops of L, B, M, A for each round
	#make a character vector of samples that want to keep
	Leps<-p_dat2[p_dat2$Species== "N. lepida",]
	Lep2<-sample(Leps$pop,2)
	Bry<-p_dat2[p_dat2$Species== "N. bryanti",]
	Bry2<-sample(Bry$pop,2)
	mac<-p_dat2[p_dat2$Species== "N. macrotis",]
	Mac2<-sample(mac$pop,2)
	Alb<-p_dat2[p_dat2$Species== "N. albigula",]
	Alb2<-sample(Alb$pop,2)

	samples<-c("Wupa_steph", "Wupa_dev", "Rio_cin", Lep2, Bry2, Mac2,Alb2)
	pp<-prune_samples(samples, ppr)

	#leaves have to be labeled identically
	df<-sample_data(pp)
	df$sp<-substr(df$code, start=1, stop=1)
	df$sp2<-make.unique(df$sp, sep = "")

	pp_dist<-phyloseq::distance(pp, "bray")
	micD<-as.matrix(pp_dist)
	row.names(micD)<-df$sp2
	colnames(micD)<-df$sp2

	micD2<-as.dist(micD)

	UP <- hclust(micD2, method="average")  #builds UPGMA tree
	micros <- as.dendrogram(UP)

	mic_new<-write.dendrogram(micros, edges = TRUE)
	mic_tree<-ape::read.tree(text=mic_new)

	distance <- TreeDistance(mic_tree, Neo_tree1)
	dist_output[i,1]<-distance
	}

distance_m<-mean(dist_output[,1])

###bootstrap p-value for 11-tipped tree pairs
perm<-100000
output<-matrix(0, nrow = perm,  ncol =1)
for(i in 1:perm){
	RT1<-RandomTree(tips=11, root = FALSE)
	RT2<-RandomTree(tips=11, root = FALSE)
	CID_RT <- TreeDistance(RT1,RT2)
	output[i,1]<-CID_RT
	}

mean(output[,1])
boots<-output[,1]
length(boots[boots> distance_m])/perm
