## Started 16 June 2025 ##
## By Avery ##

#Making a small phylogeny for Mt Rainier analysis

library(ape)
library(geiger)


#Load data
neighbordir <- "analysis/mountrainier/input/neighborhood/"
neighbordata <- read.csv(paste0(neighbordir, "neighborhooddatacored2008to2011.csv"))
neighbordata.info <- read.csv(paste0(neighbordir, "neighborhooddatacored2008to2011treeinfo.csv"))

#A little maintenance
neighbordata.info$Species <- toupper(str_trim(neighbordata.info$Species))
neighbordata$Species.neighbor. <- toupper(str_trim(neighbordata$Species.neighbor.))

neighbordata$Species.neighbor.[neighbordata$Species.neighbor. == "PSHE"] <- "PSME" #Fixing a typo
neighbordata.info$Species[neighbordata.info$Species == "PSHE"] <- "PSME" #Fixing a typo


#Phylo distance
treephylo <- read.tree("analysis/pnw/input/ALLMB.tre")
spcodes <- read.csv("analysis/mountrainier/data/sppcodes.csv")
spcodes$Species[spcodes$Species == "CANO9"] <- "CANO"
spcodes$Species[spcodes$Species == "PIMO3"] <- "PIMO"

rainiercodes <- unique(neighbordata.info$Species)
rainiercodes[!(rainiercodes %in% spcodes$Species)] #What is ACCI supposed to be? 
rainiercodes <- rainiercodes[!(rainiercodes %in% c("", "ACCI"))]

treephylo <- keep.tip(treephylo, gsub(" ","_", spcodes$Binomial[spcodes$Species %in% rainiercodes]))
write.tree(treephylo, "analysis/mountrainier/data/rainiertree.tre")