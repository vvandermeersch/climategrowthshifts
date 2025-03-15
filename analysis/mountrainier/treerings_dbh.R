## Started 11 March 2025 ##
## By Lizzie ##

## Looking for trends over time in size x growth ##

# Cribbing some off data_exploration.qmd at the start here #
# See also allometreec.R for simulation study #

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
## 
spbyplotall <- subset(treerings, select=c("Species", "Stand"))
spbyplot <- spbyplotall[!duplicated(spbyplotall),]
spbyplot <- spbyplot[with(spbyplot, order(spbyplot$Species)), ]

pdf(paste(wd, "/figures/treerings_groupbyDBH_rw.pdf", sep=""), onefile = TRUE)
for (i in c(1:nrow(spbyplot)))
	{
	speciesplothere <- treerings[which(treerings$Species==spbyplot$Species[i] & 
		treerings$Stand==spbyplot$Stand[i]),]
	plotme <- ggplot(speciesplothere, aes(x=year, y=rw, color=as.factor(Core))) +
		geom_line() +
		theme_bw() +
		xlab("year") +
		ylab("ring width (mm)") +
		ggtitle(paste(spbyplot$Species[i], "at", spbyplot$Stand[i])) +
		facet_wrap(as.numeric(DBH)~.) 
	print(plotme)	
}
dev.off()

# Now do it for basal area (aka area of a circle) and then do it for THE INCREMENT ...
# So far I just reconstrut the DBH ... need to get the two areas and take the DIFFERENCE
spbyplotplusall <- subset(treerings, select=c("Species", "Stand", "uniqueID"))
spbyplotplus <- spbyplotplusall[!duplicated(spbyplotplusall),]
spbyplotplus <- spbyplotplus[with(spbyplotplus, order(spbyplotplus$Species, spbyplotplus$Species)), ]
# Watch out! Super cheap and ugly code below!
treeringsextra <- treerings[1,]
treeringsextra$radconst <- NA
treeringsextra$rwarea <- NA
for (i in c(1:nrow(spbyplotplus)))
	{
	speciesplothere <- treerings[which(treerings$Species==spbyplotplus$Species[i] & 
		treerings$Stand==spbyplotplus$Stand[i] & 
		treerings$uniqueID==spbyplotplus$uniqueID[i]),]
	speciesplothere$radconst <- NA
	speciesplothere <- speciesplothere[with(speciesplothere, order(-speciesplothere$year)), ]
		for (j in c(1:(nrow(speciesplothere)-1)))
		{ 
			# First get the radius through time: DBH (cm to mm) divide by 2 minus increment
			speciesplothere$radconst[1] <- (as.numeric(speciesplothere$DBH[1])*100)/2
			speciesplothere$radconst[j+1] <- speciesplothere$radconst[j]-speciesplothere$rw[j]
			# Now get the area (basal area) and takes the difference over time
			speciesplothere$rwarea[1] <- NA
			speciesplothere$rwarea[j+1] <- (pi*speciesplothere$radconst[j]^2)-(pi*speciesplothere$radconst[j+1]^2)

		}	
	treeringsextra <- rbind(treeringsextra, speciesplothere)	
}
treeringsextra <- treeringsextra[-1,]
treeringsextra$radconst <- treeringsextra$radconst/100

pdf(paste(wd, "/figures/treerings_groupbyDBH_areainc.pdf", sep=""), onefile = TRUE, width=10, height=8)
for (i in c(1:nrow(spbyplot)))
	{
	speciesplothere <- treeringsextra[which(treeringsextra$Species==spbyplot$Species[i] & 
		treeringsextra$Stand==spbyplot$Stand[i]),]
	plotme <- ggplot(speciesplothere, aes(x=year, y=radconst, color=as.factor(Core))) +
		geom_line() +
		theme_bw() +
		xlab("year") +
		ylab("basal area increment (mm)") +
		ggtitle(paste(spbyplot$Species[i], "at", spbyplot$Stand[i])) + 
		facet_wrap(as.numeric(DBH)~., scales="free") 
	print(plotme)	
}
dev.off()


##
## Simulation study!
## See allometreec.R instead
## And see https://github.com/vvandermeersch/climategrowthshifts/issues/7

## This code below is working because of:
# 1) My assumption that alpha is a constant
# 2) Something bad I did with circumference (why is it negative?)
# 3) Possibly some other coding issues that I did not catch
## But I am keeping it here for now as it seemed close

# Growth goes down with increasing basal area (pi*r^2)
# Using Mike's solution for RW (ring width, aka increment)
# The below math is simplified to when: t_n-t_{n-1}=1
if(FALSE){
rwsolution <- function(circum, alpha, i){
	ringwidthhereplus <- ((-circum[i-1]/pi) + sqrt((circum[i-1]^2/pi^2)-4*alpha))/2
	ringwidthhereminus <- ((-circum[i-1]/pi) - sqrt((circum[i-1]^2/pi^2)-4*alpha))/2
	return(ringwidthhereplus)
}

# Next: Build a bunch of trees of different DBH
treen <- 20
yearn <- 100
treeid <- rep(1:treen, each=yearn)
yearz <- rep(seq(1:yearn), treen)
N <- length(yearz)
startcircum <- rnorm(treen, 30, 500)
rwdat <- rep(NA, N)
circumdat <- rep(NA, N)
alpha <- 0.002

for(n in 1:N){
	treehere <- treeid[n]
	if(yearz[n]==1){
		circumdat[n] <- startcircum[treehere]
		rwdat[n] <- 0
	} else{
		rwdat[n] <- rwsolution(circumdat, alpha, n)
		circumdat[n] <- circumdat[n-1]+rwdat[n-1]
	}
}

df <- data.frame(treeid, yearz, circumdat, rwdat)
ggplot(df, aes(x=yearz, y=rwdat, color=as.factor(treeid))) +
		geom_line() +
		theme_bw() +
		xlab("year") +
		ylab("ring width") 
ggplot(df, aes(x=yearz, y=circumdat, color=as.factor(treeid))) +
		geom_line() +
		theme_bw() +
		xlab("year") +
		ylab("circumference") 
# I have some parameter number issues ... 
# And the math is wrong?
}

