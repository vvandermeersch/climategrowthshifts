## Started 16 June 2025 ##
## By Avery ##

####Goals####
#Calculate the sum(basal area) for all trees surrounding each tree 
##-- maybe one version with all trees and one with only canopy trees
#Also calculate the sum(basal area) but weight each tree by sqrt(phylogenetic distance)
#do all of this in this repo and confirm it merges with the data that @vvandermeersch is using

library(ape)
library(geiger)
library(tidyverse)

#Load data
neighbordir <- "analysis/mountrainier/input/neighborhood/"
neighbordata <- read.csv(paste0(neighbordir, "neighborhooddatacored2008to2011.csv"))
neighbordata.info <- read.csv(paste0(neighbordir, "neighborhooddatacored2008to2011treeinfo.csv"))
rainiercodes

#A little maintenance
neighbordata.info$Species <- toupper(str_trim(neighbordata.info$Species))
neighbordata$Species.neighbor. <- toupper(str_trim(neighbordata$Species.neighbor.))

neighbordata$Species.neighbor.[neighbordata$Species.neighbor. == "PSHE"] <- "PSME" #Fixing a typo
neighbordata.info$Species[neighbordata.info$Species == "PSHE"] <- "PSME" #Fixing a typo


#Phylo distance
treephylo <- read.tree("analysis/mountrainier/data/rainiertree.tre")
spcodes <- read.csv("analysis/mountrainier/data/sppcodes.csv")
spcodes$Binomial <- gsub(" ","_",spcodes$Binomial)
spcodes$Species[spcodes$Species == "CANO9"] <- "CANO"
spcodes$Species[spcodes$Species == "PIMO3"] <- "PIMO"

rainiercodes <- unique(neighbordata.info$Species)
rainiercodes[!(rainiercodes %in% spcodes$Species)] #What is ACCI supposed to be? 
rainiercodes <- rainiercodes[!(rainiercodes %in% c("", "ACCI"))]

for(i in 1:length(treephylo$tip.label)){
  treephylo$tip.label[i] <- spcodes$Species[spcodes$Binomial == treephylo$tip.label[i]]
}

rainiercoph <- cophenetic(treephylo)

neighbordata$phylodist <- NA

#Calculate phylogenetic distance of each neighbor
for(i in 1:nrow(neighbordata)){
  mysp <- neighbordata.info$Species[neighbordata.info$Tag == neighbordata$Treeid[i] &
                                      neighbordata.info$Stand == neighbordata$Stand[i]]
  if(length(mysp)>1){neighbordata$phylodist[i] <- 999; next}
  #There are some rows where there appear to be multiple matches for both Stand and Treeid?
  neighbsp <- neighbordata$Species.neighbor.[i]
  if(length(neighbsp)>1){neighbordata$phylodist[i] <- 998; next}
  ind1 <- rownames(rainiercoph) == neighbsp
  ind2 <- colnames(rainiercoph) == mysp
  if(sum(ind1)==1 & sum(ind2)==1){
    neighbordata$phylodist[i] <- rainiercoph[ind1, ind2]
  } else{neighbordata$phylodist[i] <- NA}
}

#Set up data frame
area.df <- as.data.frame(unique(neighbordata[,c("Stand","Treeid")]))
area.df$neighborarea <- NA
area.df$neighborarea.canopy <- NA
area.df$phylowarning <- NA

neighbordata$phylodist.trans <- NA
neighbordata$phylodist.trans[neighbordata$phylodist==0] <- 1
neighbordata$phylodist.trans[neighbordata$phylodist!=0 & !is.na(neighbordata$phylodist)] <-
  1/sqrt(neighbordata$phylodist[neighbordata$phylodist!=0 & !is.na(neighbordata$phylodist)])

#Treeids are not replicated across stands, so just do analysis on Treeid
for(i in 1:nrow(area.df)){
  tmp <- neighbordata[neighbordata$Treeid == area.df$Treeid[i],]
  area.df$neighborarea[i] <- sum(tmp$DBH.neighbor., na.rm = T)
  area.df$neighborarea.canopy[i]<- sum(tmp$DBH.neighbor.[tmp$Canopy. == "y"])
  
  area.df$neighborarea.phylo[i] <- sum(tmp$DBH.neighbor. * tmp$phylodist.trans)
  area.df$neighborarea.canopy.phylo[i] <- sum(tmp$DBH.neighbor.[tmp$Canopy. == "y"] *
                                                tmp$phylodist.trans)
  if(any(tmp$phylodist >= 998 | is.na(tmp$phylodist))){
    area.df$phylowarning[i] <- T
  } #Tagging things with phylogenetic distances messed up due to non-singular tag/stand matching
}

hist(log(area.df$neighborarea), xlab ="log Neighbor Basal Area")
hist(log(area.df$neighborarea.phylo), xlab ="log Neighbor Basal Area (Canopy)")
hist(log(area.df$neighborarea.canopy), xlab ="log Neighbor Basal Area (Phylo-weighted)")
hist(log(area.df$neighborarea.canopy.phylo), xlab ="log Neighbor Basal Area (Canopy, Phylo-weighted")
