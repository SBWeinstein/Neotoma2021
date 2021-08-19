# script for stable isotope mixing models in R 
# R version 3.5.3 (2019-03-11)
library(ggplot2)
library(simmr) # v.0.4.2

# will be using δ13C reference values for C3, C4, and CAM plants from:
# O'Leary, M. H. 1988. Bioscience: 38

# set working directory:
setwd("<directory_location>")

# input files
SIdata <- read.csv("SIphylodata.csv") # woodrat fecal δ13C data
IDs <- read.csv("isotopeIDs.csv") # IDs for each sample
s_means = read.csv("C3-reference-values.csv", header=FALSE, sep=",") # isotope reference values
s_sds = read.csv("C3-reference-valuesSD.csv", header=FALSE, sep=",") # isotope reference SD values

# formatting woodrat fecal δ13C data
mix = as.matrix(SIdata)

# formatting reference plant values (C3, CAM, C4)
colnames(mix) = c('δ13C')
s_names=c('C3','CAM','C4')
s_means = as.matrix(s_means)
s_sds = as.matrix(s_sds)

# Loading the data matrix into simmr:
simmr_in = simmr_load(mixtures=mix, source_names = s_names, source_means = s_means, source_sds = s_sds)

# tracer plot
plot(simmr_in)

# running the mixing model to estimate dietary proportions
simmr1 = simmr_mcmc(simmr_in)

# Analyzing the model: confidence intervals, SDS, and quartile information
summary(simmr1,type='diagnostics')
summary(simmr1,type='statistics')
summary(simmr1,type='quantiles')

# check the fit of the model with a posterior regression
posterior_predictive(simmr1)

# look at graphical summaries of the model run
prior_viz(simmr1) # shows the prior vs posterior predictions

# Plot the density of the isotopic signatures against their relative proportions compared to each reference
plot(simmr1,type='density')

# Histogram-Contour plots:
# This shows the source histograms on the diagonal, contour plots of the relationship between the sources on
# the upper diagonal, and the correlation between the sources on the lower diagonal.
plot(simmr1,type='matrix') # this shows that C4 vs CAM signatures overlap and the model cannot disintquish between these two sources

# comparing theoretical dietary sources for the experimental samples
compare_sources(simmr1)

# Combine C4 and CAM sources as the model cannot discriminate between the two:
simmr_combine = combine_sources(simmr1,to_combine=c('CAM','C4'),new_source_name='CAM/C4')

# create a new tracer plot with the combined sources:
plot(simmr_combine$input)

# combined density plot
plot(simmr_combine,type='density')

# re-run the model to compare sources with the combined Cam/C4 dietary sources
compare_sources(simmr_combine)

# To generate individual source proportions for each animal:
simmr_groups = simmr_load(mixtures=mix,
                          source_names=s_names,
                          source_means=s_means,
                          source_sds=s_sds,
                          group=as.factor(paste('ID', IDs)))

simmr_groups_out = simmr_mcmc(simmr_groups,mcmc_control=list(iter=10000,burn=1000,thin=10,n.chain = 4))

simmr_groups_combined = combine_sources(simmr_groups_out,to_combine=c('CAM','C4'),new_source_name='CAM/C4')
compare_sources(simmr_groups_combined)

# write out an summaries of dietary proportions and SD for each individual sample
out <- summary(simmr_groups_combined,type='statistics',group=c(1:63))
out2 <- (out[["statistics"]])
write.csv(out2, "phylo-isotope-proportions.csv")
