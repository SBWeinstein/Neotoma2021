# Neotoma2021 
## This respository contains R code associated with "Microbiome stability and structure is governed by host phylogeny over diet and geography in woodrats (*Neotoma* spp.)" 

### R scripts for manuscript are as follows:

**16s_processing.R:** Initial processing of bacterial 16s amplicon sequences,
intitial filtering, and quality control. Inputs are fastq.gz files from each sequencing run.
outputs are rarefied and non-rarefied microbiome phyloseq objects for downstream analyses.

**trnl_processing.R:** Initial processing of plant trnL amplicon sequences.
This script is very similar to the 16s_processing.R file, with modifications
for different primers and different taxonomy assignment methods.
Taxonomy is assigned using separate custom python script: see https://github.com/robertgreenhalgh/stand. 
Inputs are fastq.gz files from each sequencing run. Output is a plant ASV table. 

**cactus_adjust.R:** Adjust sequencing based diet data to account for missing cactus reads using stable isotope data. 
Script creates a phyloseq object from trnl data, filters ASVs, and creates basic plots demonstrating missing cactus in diets. 
Individual diets are adjusted using carbon stable isotope data based on SIMMR models from captive cactus diet trials and wild collected fecal samples. Produces the cactus adjusted, filtered, rarefied plant OTU table used in downstream analyses.

**Diversity_Composition_comps.R:** Script for alpha and beta diversity analyses for wild rat microbiome and diet data. 
Includes code for Figure 2.

**MRM_models.R:** Script for multiple regression on distance matrices (MRM) models and code for Figure 3.
Code for Figure 4B is also based on this script.

**Differential_Abundance_diet.R:** Script for differential abundance analyses identifying ASVs (or bacterial families)
associated with different wild diet components using DESeq2.

**Sympatry.R:** Script for testing whether heterospecifics from the same site had more similar
microbiomes (based on Jaccard similarity) than matched species pairs from different
sites. Script also includes code to compare average percent of shared ASVs for sympatric and
equivalent allopatric populations, and code for Figure S5.

**mantel_and_clustering.R:** Code for two approaches for examining correlations between
the microbiome and geography, phylogeny, or diet. Includes comparisons at the individual,
population, and species levels.  Mantel and Partial Mantel tests are based on distance
matrices. Hierarchical clustering analyses compare tree congruence, calculating significant
congruence via a bootstrapped p-value.

**Differential_Abundance_captivity.R:** Code for differential abundance analyses identifying ASVs (genera, and families)
that increased/decrease in wild v. captive animals. Includes analyses of taxa that change in
each population and when all individual are grouped together.
Script is very similar to Differential_Abundance_diet.R but includes loop to run
through all populations.

**Captivity_impacts.R:** Exploring how captivity impacts microbiome diversity and
composition. Includes PERMANOVAs, builds a dataframe of changes in observed richness and composition for each animal, 
and includes statistical analyses using that dataframe.

**Morans.R:** Testing for phylogenetic signal in how animal microbiomes change in captivity.
Code uses the dataframe produced in Captivity_impacts.R.

**rare_lost.R:** Using paired wild and captive 16s data to test whether ASVs that
are less abundant in wild animals are more likely to be lost in captivity.

**Captivity_homogenize.R:** Code to test whether captivity reduces variation
among individuals.

**PNM_sites_taxa.R:** Code for fitting the prokaryote neutral model (PNM) to ASV data,
testing whether model fit correlates with host species diversity at a site,
calculating model fit for wild and captive animals, and identifying "neutral"
and "non-neutral" taxa. The code for the PNM comes from the supplemental material (Supplemental Code 1) from
Burns, A., Stephens, W., Stagaman, K. et al. Contribution of neutral processes
to the assembly of gut microbial communities in the zebrafish over host development.
ISME J 10, 655–664 (2016). https://doi.org/10.1038/ismej.2015.142.
