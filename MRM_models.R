#distance matrix regression,
#includes loop to output values for BC, J, WU, U distances from  microbiome data

setwd("<>") #setworking directory

library(ggplot2)
library(phyloseq)
library(phangorn)
library(tidyr)
library(geosphere)
library(ecodist)
library(eulerr)#for euler diagrams in figures
library(dplyr) #use to round values in MRM output table

#files
Neo_tree<-read.tree(file="Neotoma_phylosymtree.nwk") #Tree from Matoqc
site<-as.data.frame(read.csv("Site_lat_long.csv")) #lat/long data
plant<-readRDS( file = "diet1percent_rar_29May20.rds") #diet data
ps.rout<-readRDS("ps.rout_phylo_18Oct2020.rds") #microbiome data

#use phylogeny, diet, geography data to make distance matrices for MRM models
#make distance matrix for phylogeny using rat tree
#replace tip labels to match format of species names in microbiome data
Neo_tree$tip.label<-c("N. albigula", "N. stephensi", "N. macrotis",
                      "N. cinerea", "N. bryanti", "N. lepida", "N. devia" )
#compute the pairwise distances between the pairs of tips using branch lengths
D_evo<-cophenetic(Neo_tree)
#convert into long format, species1, species2, distance
D_evo_long<-gather(as.data.frame(D_evo), key="Species2", value= "distance")
D_evo_long["Species1"]<-c(rep(rownames(D_evo), 7)) #make extra column for compared to species
D_evo_long <- D_evo_long[c(3,1,2)]  #reorder columns

#make distance matrix for geography using lat/long data
D_site<-distm(cbind(site$longitude,site$latitude), fun=distGeo) #longitude first
row.names(D_site) <- site$site
colnames(D_site) <- site$site
D_site_km<-D_site/1000 #convert to km
#convert to long format
D_site_long<-gather(as.data.frame(D_site_km), key="Site2", value= "distance")
D_site_long["Site1"]<-c(rep(rownames(D_site_km), 19))
D_site_long <- D_site_long[c(3,1,2)]

#Make  diet distance matrix
#using individual diet data, at family level
psF<-tax_glom(plant, "Family")
D_bray_diet<-phyloseq::distance(psF, "bray")
dietm<-as.matrix(D_bray_diet)
row.names(dietm) <- sample_data(psF)$Woodrat.ID
colnames(dietm) <- sample_data(psF)$Woodrat.ID
#convert to long format
D_diet_long<-gather(as.data.frame(dietm), key="pop2", value= "distance")
D_diet_long["pop1"]<-c(rep(rownames(dietm), dim(dietm)[1]))
D_diet_long <- D_diet_long[c(3,1,2)] #reorder columns

#######################
#create microbiome distance matrix
#subset to wild (repeat with subset to captive)
psW<-subset_samples(ps.rout, Sample_type=="wild") #or "captive"
psW<-prune_taxa(taxa_sums(psW) > 0, psW) #2543 in all captives, 2704 wild

#restrict to animals with diet data
psW<-subset_samples(psW, Woodrat.ID %in% row.names(dietm))
psW<-prune_taxa(taxa_sums(psW) > 0, psW)# wild 2697 taxa and 147 samples

##################################
#for neutral/non-neutral taxa MRM output comparisons,
#subset taxa table based on neutral model outputs, and run models with subsets
#for example: for the list of "neutral" taxa
#psN<-subset_taxa(psAll, rownames(tax_table(psAll)) %in% neutral )
####################################

#MRM models
ps_data<-sample_data(psW)

###Make big host phylogeny matrix to match microbiome dimensions
#start with blank new matrix with dims and names matching the microbiome matrix
my.matrix <- matrix(0, nrow = length(ps_data$Species), ncol = length(ps_data$Species))
row.names(my.matrix) <- ps_data$Species
colnames(my.matrix) <- ps_data$Species

#use the table to fill in values in matrix with dimensions matching all the samples
for (i in 1:nrow(D_evo_long)) {
   my.matrix[rownames(my.matrix) == D_evo_long$Species2[i],
	 colnames(my.matrix) == D_evo_long$Species1[i]] <- D_evo_long$distance[i]
}
big_evo_matrix<-my.matrix #rename

#Make big SITE matrix to match microbiome dimensions
#make blank new matrix with dims and names matching the microbiome matrix
big.site.matrix <- matrix(0, nrow = length(ps_data$site), ncol = length(ps_data$site))
row.names(big.site.matrix) <- ps_data$site
colnames(big.site.matrix) <- ps_data$site

#use the table to fill in values in matrix with dimensions matching all the samples
for (i in 1:nrow(D_site_long)) {
   big.site.matrix[rownames(big.site.matrix) == D_site_long$Site2[i],
	 colnames(big.site.matrix) == D_site_long$Site1[i]] <- D_site_long$distance[i]
}

#make big DIET matrix with dims and names matching the microbiome matrix
big.diet.matrix <- matrix(0, nrow = length(ps_data$Woodrat.ID), ncol = length(ps_data$Woodrat.ID))
row.names(big.diet.matrix) <- ps_data$Woodrat.ID
colnames(big.diet.matrix) <- ps_data$Woodrat.ID

#use the table to fill in values in matrix with dimensions matching all the samples
for (i in 1:nrow(D_diet_long)) {
   big.diet.matrix[rownames(big.diet.matrix) == D_diet_long$pop2[i],
	 colnames(big.diet.matrix) == D_diet_long$pop1[i]] <- D_diet_long$distance[i]
}

###################################################
###loop to calculate values for all distance metrics
dist_list<-c("bray", "jaccard", "wunifrac",  "unifrac")
output<-matrix(0, nrow =length(dist_list),  ncol =19 )
colnames(output)<-c("dist", "DPSr", "DPr", "DSr", "PSr", "Sr", "Dr", "Pr",
					"Diet","Phy", "Site",
				"PhyDiSi", "DietPhy", "PhySite", "DietSite","DPSp","Sp", "Dp", "Pp")

for(j in 1:length(dist_list)){
	output[j,1]<-dist_list[j]
	D_bray_ps<-phyloseq::distance(psW, method=dist_list[j])
	dm<-as.matrix(D_bray_ps)

	#matrices need to be a "distance matrix"
	DPS<-MRM(as.dist(dm) ~ as.dist(log10(big.site.matrix+1)) +
          as.dist(big_evo_matrix)+ as.dist(big.diet.matrix),  nperm=1000)

	#For variance partioning
	DP<-MRM(as.dist(dm) ~  as.dist(big_evo_matrix)+ as.dist(big.diet.matrix),  nperm=1000)
	DS<-MRM(as.dist(dm) ~  as.dist(log10(big.site.matrix+1))+ as.dist(big.diet.matrix),  nperm=1000)
	PS<-MRM(as.dist(dm) ~  as.dist(log10(big.site.matrix+1))+ as.dist(big_evo_matrix),  nperm=1000)
	S<-MRM(as.dist(dm) ~ as.dist(log10(big.site.matrix+1)),  nperm=1000)
	D<-MRM(as.dist(dm) ~ as.dist(big.diet.matrix),  nperm=1000)
	P<-MRM(as.dist(dm) ~ as.dist(big_evo_matrix), nperm=1000)

	DPSr<-DPS$r.squared[1] #r values are 1
	DPSp<-DPS$r.squared[2] #p values are 2
	DPr<-DP$r.squared[1]
	DSr<-DS$r.squared[1]
	PSr<-PS$r.squared[1]
	Sr<-S$r.squared[1]
	Sp<-S$r.squared[2]
	Dr<-D$r.squared[1]
	Dp<-D$r.squared[2]
	Pr<-P$r.squared[1]
	Pp<-P$r.squared[2]

	output[j,2]<-DPS$r.squared[1]
	output[j,3]<-DP$r.squared[1]
	output[j,4]<-DS$r.squared[1]
	output[j,5]<-PS$r.squared[1]
	output[j,6]<-S$r.squared[1]
	output[j,7]<-D$r.squared[1]
	output[j,8]<-P$r.squared[1]

	#get contributions
	output[j,9]<-DPSr-PSr
	output[j,10]<-DPSr-DSr
	output[j,11]<-DPSr-DPr
	#overlaps
	output[j,12]<- DPSr +Dr + Pr+ Sr-DPr-DSr-PSr
	output[j,13]<- PSr +  DSr - Sr - DPSr
	output[j,14]<- DSr +DPr - Dr -DPSr
	output[j,15]<- PSr + DPr - Pr-DPSr
	#pvalues
	output[j,16]<- DPSp
	output[j,17]<- Sp
	output[j,18]<- Dp
	output[j,19]<- Pp
	}

op1<-data.frame(output)

#library(dplyr)  # to round values in table
op2<-mutate_if(op1[,-1], is.factor, ~ as.numeric(as.character(.x)))
op3<-round(op2,digits=3)
op4<-data.frame(metric=op1[,1], op3)

########Make slope graphs for figure 3, make each graph
op_wild<-op4
op_wild$type<-rep("W", 4)

#repeat for captive animal data
op_cap<-op4
op_cap$type<-rep("C", 4)
vp<-rbind(op_cap, op_wild) #make single dataframe with wild and captive outputs

#subset for plots with each factor
vp1<-vp[,c(1,2,6,7,8, 20)]
ALL<-vp1[,c(1,2,6)]
SITE<-vp1[,c(1,3,6)]
DIET<-vp1[,c(1,4,6)]
PHY<-vp1[,c(1,5,6)]

DATA<-ALL #pick factor
colnames(DATA)<-c("metric", "var", "type")
colpick<-"black"	#"black" "darkgreen","darkmagenta", "blue3"

plot<-ggplot(data = DATA, aes(x = reorder(type, desc(type)), y = var, group = metric)) +
	geom_line(aes(linetype = metric), size = 2, color=colpick, alpha=.8) +
  geom_point( size = 4, color=colpick, alpha=.8) +
  ylim(0, .63)+
 	theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
 	theme(axis.text.x = element_text(size=30)) +
	theme(axis.text.y = element_text(size=18)) +
  theme(axis.ticks.x = element_blank())

############################
##to make Euler diagrams. Note From Legendre 2008:
#Negative values of partial R2 are interpreted as zeros;they correspond to cases
#where the explanatory variables explain less variation than random normal variables would.

#make wild venn w bray curtis, repeat for captive
Bray<-op_wild[1,]
#for making venn diagram
D100<-as.vector(Bray$Diet)*100
P100<-as.vector(Bray$Phy)*100
S100<-as.vector(Bray$Site)*100
PS100<-as.vector(Bray$PhySite)*100
PD100<-as.vector(Bray$DietPhy)*100
DS100<-as.vector(Bray$DietSite)*100
PDS100<-as.vector(Bray$PhyDiSi)*100

##plot as venn diagram
#the R circle gives a reference size for 100% of variance
#could do math and make it overlap now. but math is hard.  center in Illustrator
combo_wild <- c(Diet = D100, Phylogeny = P100, Site = S100,
	"Diet&Phylogeny" = PD100, "Diet&Site" = DS100, "Phylogeny&Site" = PS100,
	"Diet&Phylogeny&Site"=PDS100, "R"=100, "R&Phylogeny"=.0001)

ven<-plot(euler(combo_wild, input = c("disjoint"),shape =  "ellipse"),
	fills = list(fill = c("darkgreen","darkmagenta", "blue3" ), alpha = 0.8),
	quantities = FALSE)
