## Started 11 March 2025 ##
## By Lizzie ##

## Cribbing some off data_exploration.qmd ##

# housekeeping 
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# libraries
library(ggplot2)

wd <- "~/Documents/git/projects/grephon/climategrowthshifts/analysis/mountrainier"

treerings <- read.csv(file.path(wd, "data", "treerings_ailene", "SouthSidecoresX.csv"), header=TRUE)
treerings$uniqueID <- paste(treerings$Tag, treerings$Core, sep = "_")                      
meas_cols <-  names(treerings)[grep("X", names(treerings))]
treerings <- reshape(treerings, dir = "long", 
                     varying = meas_cols, 
                     times = as.numeric(gsub("X", "", meas_cols)),
                     v.names = "rw", timevar = "year")
treerings <- treerings[c("Stand", "Species", "Tag", "Core", "uniqueID", "DBH", "RingCount", "year", "rw")]
treerings <- na.omit(treerings)
treerings <- treerings[order(treerings$year),]
treerings$cumrw <- ave(treerings$rw, treerings$uniqueID, FUN = cumsum)
unique(treerings$Stand)

## Okay, let's figure out some size class groups to look at...
dbhspp <- subset(treerings, select=c("Stand", "Species", "DBH"))
dbhspp <- dbhspp[!duplicated(dbhspp),]
dbhspp$DBH <- as.numeric(dbhspp$DBH)

ggplot(dbhspp, aes(x=DBH, group=Species)) +
	geom_histogram() +
	facet_grid(Species~.)

# ABAM or TSHE would be good for shade tolerant
ggplot(subset(dbhspp, Species=="Abam"), aes(x=DBH, group=Stand)) +
	geom_histogram() +
	facet_grid(Stand~.)

# And for shade-intolerant
ggplot(subset(dbhspp, Species=="Psme"), aes(x=DBH, group=Stand)) +
	geom_histogram() +
	facet_grid(Stand~.)

# AE10 has a ton of size classes, so let's start there for Abam
# For Tshe: TO04 but order by DBH?

abamae10 <- subset(treerings, Species=="Abam" & Stand=="AE10")
ggplot(abamae10, aes(x=year, y=rw, color=as.factor(Core))) +
	geom_line() +
	theme_bw() +
	ggtitle("ABAM at AE10") +
	facet_wrap(as.numeric(DBH)~.) 

quartz()
tsheto04 <- subset(treerings, Species=="Tshe" & Stand=="TO04")
ggplot(tsheto04, aes(x=year, y=rw, color=as.factor(Core))) +
	geom_line() +
	theme_bw() +
	ggtitle("TSHE at TO04") +
	facet_wrap(as.numeric(DBH)~.) 

quartz()
psmeto04 <- subset(treerings, Species=="Psme" & Stand=="TO04")
ggplot(psmeto04, aes(x=year, y=rw, color=as.factor(Core))) +
	geom_line() +
	theme_bw() +
	ggtitle("PSME at TO04") +
	facet_wrap(as.numeric(DBH)~.) 

## Next up, just build a loop for each species (PDF) x stand...



	
