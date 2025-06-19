## Started 10 June 2025 ##
## By Lizzie but following grephon code, which follow coes I got from Deirdre ##

# housekeeping
rm(list=ls()) 
options(stringsAsFactors=FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
if(length(grep("lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/grephon/climategrowthshifts/analysis/pnw")
} else if(length(grep("victor", getwd())>0)) { 
  setwd("/home/victor/projects/climategrowthshifts/analysis/pnw")
}  else if (length(grep("helloworld", getwd()))>0) {
  setwd("boomboom")
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

# Make the tree ultrametric
run <- FALSE
if(run){
  
  library(chronos)
  library(doFuture)
  
  .droptip <- function(phy, tip) {
    s <- which(phy$edge[, 2] == tip)
    phy$edge.length <- phy$edge.length[-s]
    phy$edge <- phy$edge[-s, ]
    k <- which(phy$edge > tip)
    phy$edge[k] <- phy$edge[k] - 1L
    phy$tip.label <- phy$tip.label[-tip]
    phy
  }
  
  CVmod <- function(phy, lambda = 10^(-3:2), model = "correlated",
                    quiet = FALSE, calibration = makeChronosCalib(phy))
  {
    cv <- numeric()
    n <- Ntip(phy)
    nb <- n - 1L
    D <- numeric(n)
    one2n <- 1:n
    S <- match(one2n, phy$edge[, 2])
    cal2 <- calibration
    cal2$node <- cal2$node - 1L
    
    ctrl <- chronos.control(dual.iter.max = 400) # added by VV
    
    for (l in lambda) {
      if (!quiet)
        cat("Running cross-validation with lambda = ", l, "\n", sep = "")
      chr <- chronos(phy, l, model, quiet = TRUE, calibration = calibration, control = ctrl)
      for (i in one2n) {
        if (!quiet)
          cat("\r    Dropping terminal branch ", i, "/", n, sep = "")
        trb <- .droptip(phy, i)
        chrb <- chronos(trb, l, model, quiet = TRUE, calibration = cal2, control = ctrl)
        s <- S[i]
        ancb <- phy$edge[s, 1L] - 1L
        ## ancb is the node# of the direct ancestor of the dropped tip *in* the modified tree
        btb <- branching.times(chrb)
        j <- if (ancb == nb + 1L) 1L else 2L
        xstar <- btb[ancb - nb] * attr(chrb, "rates")[which(trb$edge[, j] == ancb)]
        D[i] <- (xstar - phy$edge.length[s])^2 / xstar
      }
      sD <- sum(D)
      if (!quiet) cat("   ==>   CV = ", sD, "\n", sep = "")
      cv <- c(cv, sD)
    }
    cbind(lambda = lambda, CV = cv)
  }
  
  # The idea is to find the lambda values with the lowest CV
  b <- Sys.time()
  plan(multisession, workers = 2)
  lambdas <- 10^seq(-6, 0, 0.5)
  res <- foreach(lambda = lambdas, .combine=rbind, .options.future = list(seed = TRUE)) %dofuture% {
    
    resl <- CVmod(phy.plants.here, lambda)
    resl
  }
  plan(sequential);gc()
  e <- Sys.time()
  plot(res, type = "o", log = "xy")
  
  # We can narrow the search
  b <- Sys.time()
  plan(multisession, workers = 4)
  lambdas <- 10^seq(-6, -5, 0.1)
  res <- foreach(lambda = lambdas, .combine=rbind, .options.future = list(seed = TRUE)) %dofuture% {
    
    resl <- CVmod(phy.plants.here, lambda)
    resl
  }
  plan(sequential);gc()
  e <- Sys.time()
  print(e-b)
  plot(res, type = "o", log = "xy")
  
  # And narrow again?
  b <- Sys.time()
  plan(multisession, workers = 4)
  lambdas <- 10^seq(-5.6, -5.4, 0.02)
  res <- foreach(lambda = lambdas, .combine=rbind, .options.future = list(seed = TRUE)) %dofuture% {
    
    resl <- CVmod(phy.plants.here, lambda)
    resl
  }
  plan(sequential);gc()
  e <- Sys.time()
  print(e-b)
  plot(res, type = "o", log = "xy")
  
  # And narrow again?
  b <- Sys.time()
  plan(multisession, workers = 3)
  lambdas <- 10^seq(-5.52, -5.48, 0.005)
  res <- foreach(lambda = lambdas, .combine=rbind, .options.future = list(seed = TRUE)) %dofuture% {
    
    resl <- CVmod(phy.plants.here, lambda)
    resl
  }
  plan(sequential);gc()
  e <- Sys.time()
  print(e-b)
  plot(res, type = "o", log = "xy")
  
  
  #Once you choose lambda, you can make the tree ultrametric with chronos
  
  # Make it ultrametric
  ctrl <- chronos.control(dual.iter.max = 1e6)
  phy.plants.here.ultra <- chronos(phy.plants.here, lambda =  3.162278e-06, model = "correlated", control = ctrl)
  #log-Lik = -46.6315
  #PHIIC = 157.26
  
  #Resolve multichotomies
  phy.plants.here.ultra <- multi2di(phy.plants.here.ultra)
  
  # Rescale with the fossil information of how old is the phylogeny (My Botryo family was 69.9 Mya)
  # treeFam_ultraRes <- rescale(treeFam_ultraR, model = "depth", 69.9)
  
}
