## Started 16 June 2025 ##
## By Avery ##

####Goals####
#Calculate the sum(basal area) for all trees surrounding each tree 
##-- maybe one version with all trees and one with only canopy trees
#Also calculate the sum(basal area) but weight each tree by sqrt(phylogenetic distance)
#do all of this in this repo and confirm it merges with the data that @vvandermeersch is using

#Load data
neighbordir <- "analysis/mountrainier/input/neighborhood/"
neighbordata <- read.csv(paste0(neighbordir, "neighborhooddatacored2008to2011.csv"))
neighbordata.info <- read.csv(paste0(neighbordir, "neighborhooddatacored2008to2011treeinfo.csv"))

#Set up data frame
area.df <- as.data.frame(unique(neighbordata[,c("Stand","Treeid")]))
area.df$neighborarea <- NA
area.df$neighborarea.canopy <- NA

#Treeids are not replicated across stands, so just do analysis on Treeid
for(i in 1:nrow(area.df)){
  tmp <- neighbordata[neighbordata$Treeid == area.df$Treeid[i],]
  area.df$neighborarea[i] <- sum(tmp$DBH.neighbor., na.rm = T)
  area.df$neighborarea.canopy[i]<- sum(tmp$DBH.neighbor.[tmp$Canopy. == "y"])
}
