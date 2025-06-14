## Started 10 June 2025 ##
## By Lizzie but following grephon code, which follow coes I got from Deirdre ##

# housekeeping
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/grephon/climategrowthshifts/analysis/pnw")
} else if (length(grep("helloworld", getwd()))>0) 
{setwd("boomboom")
}

library(stringr)
library(ape)
library(phytools)
library(geiger)

colnames <- c("latbi", "shortname")
sppfull <- read.csv("input/sppnames.csv") # I literally copied this out of https://github.com/vvandermeersch/climategrowthshifts/issues/25
names(sppfull) <- colnames

temp <- str_split_fixed(sppfull$latbi, " ", 3)
sppfull$phylo.name <- paste(temp[,1], temp[,2], sep="_")

spp <- sppfull
spp$genus <- temp[,1]
spp$species <- temp[,2]

sps.list <- sort(unique(spp$phylo.name))
genus.list=sort(unique(spp$genus))

## load phylo (from Smith and Brown 2019)
phy.plants <- read.tree("input/ALLMB.tre")
is.ultrametric(phy.plants) # we need to fix this, I am waiting on some code

## getting a list of genera in S&B's phylo
phy.genera <- unlist(
  lapply(strsplit(phy.plants$tip.label, "_"),function(x){return(x[1])}))

phy.genera.uniq<-sort(unique(phy.genera))

## how many grephon species are in the phylogeny?
phenosp.genus.inphylo <- genus.list[which(genus.list %in% phy.genera.uniq)]

# All genera are in tree
length(phenosp.genus.inphylo)
length(unique(spp$genus))

## first prune the phylogeny to include only these genera
phy.genera.here <- drop.tip(phy.plants,
                            which(!phy.genera %in% phenosp.genus.inphylo)) 

rm(phy.plants)
View(sort(phy.genera.here$tip.label))

# What's not going to merge? Nothing. Wow. 
# So I am deleting my code that deals with this, but can add it later
# See plotspeciesphylo.R in grephon repo. 
sps.list[which(!sps.list %in% phy.genera.here$tip.label)]

# now prune just the species I want
phy.plants.here <-  drop.tip(phy.genera.here,
                            which(!phy.genera.here$tip.label %in% sps.list))

# Some QUICK plots ... 
plot(phy.plants.here,cex=1.25)
plot(phy.plants.here,cex=1.25, type="f")

