#Code for simulating tree rings based on neighbor competition and phylogenetic distance
#And runs a stan model, which seems effective
#Created by Avery, June 25, 2025

library(ape)
library(geiger)
library(phytools)
library(rstan)
library(scales)

source("mcmc_analysis_tools_rstan.R")
source("mcmc_visualization_tools.R")


treephylo <- read.tree("analysis/mountrainier/data/rainiertree.tre")
rainiercoph <- cophenetic(treephylo)

sp.nums <- 1:length(treephylo$tip.label)
names(sp.nums) <- treephylo$tip.label

#Params----
#Setting parameters for simulation
Beta1 <- runif(1, 0, 5) #Beta acting on everything
alpha1 <- 0.0005#runif(1, 0, 0.1) #Exponent for BAnb
BA0 <-10 #Baseline basal area for BAnb
kappa <- c(0.95, 0.999) #Bounds for phylogenetic effect

#Base growth params
BAfmdpt <- 500
BAfmin <- 3
BAfvar <- 400
BAfB <- 4
BAfmax <- 10
BAf.kappa <- 10^-2

#Functions----
#Base growth function from basal area of focal tree
BAf.grow <- function(BAf.var, 
                     BAfmax.var = BAfmax, min=BAfmin,
                     midpt = BAfmdpt, k = BAf.kappa){
  #exp(-(BAf.var-midpt)^2/growthvar)*BAfB.var+min
  
  #Logistic growth
  (BAfmax.var-BAfmin)/(1+exp(-k*(midpt-BAf.var))) + BAfmin
}

#Cost imposed by difference in traits with neighbor, as ratio to original growth
Trait.cost <- function(td, k=kappa, beta = 1){
  1-k*exp(-beta*abs(td))}

#Outcome of competition on growth: Trait cost, multipled by neighbor Basal area function
Comp.fun <- function(BAnb, TDnb, 
                     alpha = 0.99, BA0 = 10, 
                     beta.phylo = 1, kappa = 0.05){
  lnb <- length(BAnb)
  prod(unlist(sapply(1:lnb, function(x){
    ((BAnb[x]/BA0)^-alpha)*Trait.cost(TDnb[x], beta = beta.phylo, kappa)}
  )))
}

#Simulation----
#Simulating phylo effect on basal area
nsim <- 1000
tmpneighb <- neighbordata[!is.na(neighbordata$phylodist),]
tmpneighb$BAnb <- (tmpneighb$DBH.neighbor./2)^2
tmpneighb$sp.nums <-  sp.nums[tmpneighb$Species.neighbor.]
tmpneighb <- tmpneighb[,-which(colSums(is.na(tmpneighb))==nrow(tmpneighb))]
tmpneighb <- tmpneighb[-which(rowSums(is.na(tmpneighb))>0),]

#Simulate base areas of focal species
#Choose random species, sample DBH from observed
spsim <- sample(treephylo$tip.label, size = nsim, replace = T)
DBHf<- sample(neighbordata.info$DBH[!is.na(neighbordata.info$DBH)], size = nsim, replace = T)
BAf <- (DBHf/2)^2

Nneighb <- sample(sapply(unique(neighbordata$Treeid), function(x){sum(neighbordata$Treeid==x)}),
                  nsim, replace = T)

#Simulate latent traits
latenttraits <- fastBM(treephylo)

#Simulate neighbors
Neighbors <- list()
for(i in 1:nsim){
  Neighbors[[i]] <- tmpneighb[sample(1:nrow(tmpneighb), Nneighb[i]),]
  for(j in 1:nrow(Neighbors[[i]])){
    Neighbors[[i]]$phylodist[j] <- latenttraits[spsim[i]]-latenttraits[Neighbors[[i]]$Species.neighbor.[j]]
  }
  if(i %% 1000 == 0)(message(i, "..."))
}

#Simulate data
mu = rnorm(max(sp.nums), 0.6, 0.1)
base_var = runif(max(sp.nums), 0.2, 0.4)
alpha1 = 0.01
kappa <- runif(1, 0.01, 0.07)
Beta1 <- runif(1, 0.8, 1.2)
Beta2 <- runif(1, 0.8, 1.2)
BA0 <- 10
treesp <- as.integer(sp.nums[spsim])
BaseGrowth <- unlist(sapply(1:nsim, function(x){rlnorm(1, 
                                                mu[treesp[x]],
                                                base_var[treesp[x]])}))
AllComps <- unlist(sapply(1:nsim, function(x){
  Comp.fun(Neighbors[[x]]$BAnb, Neighbors[[x]]$phylodist,
           alpha = alpha1, beta.phylo = Beta2, kappa = kappa,
           BA0 = BA0)
}))
RingGrowth <- Beta1*BaseGrowth*AllComps
hist(RingGrowth, breaks = seq(0,ceiling(max(RingGrowth)), 0.25))

#Demo rings
demorings <- read.csv("analysis/mountrainier/data/treerings_ailene/SouthSidecoresX.csv")
#hist(log(demorings$X2004))

#Basic Stan model
library(rstan)

#Ok, theres a log normal thing
#Each observation is an individual of a species, a ring width, and 
#for each individual, there are some neighbors
#Each neighbor has a size, and is a species some phylo distnace from focal species

#Lognormal dist has a mean and variance
#There is an alpha (for BA neighbor)
#A BA0 (not a parameter?)
#And kappa1 and kappa2, scaling phylogeny

start_ids <- c(1, 1+cumsum(Nneighb[-length(Nneighb)]))
end_ids <- cumsum(Nneighb)

standata <- list(
  N_trees = nsim,
  S = max(sp.nums),
  y = RingGrowth,
  N_neighbors = sum(Nneighb),
  tree_sp = treesp,
  tree_N_neighbs = Nneighb,
  tree_start_idxs = start_ids,
  tree_end_idxs = end_ids,
  
  neighbor_BA = unlist(c(sapply(1:nsim, function(x){Neighbors[[x]]$BAnb}))),
  neighbor_sp = unlist(c(sapply(1:nsim, function(x){Neighbors[[x]]$sp.nums}))),
  BA0 = 10,
  Sigma = vcv(treephylo)
  )

options(mc.cores = parallel::detectCores())
testfit <- stan("analysis/mountrainier/stan/avery_phylo.stan",
     model_name = "testrun!",
     data = standata, 
     chains = 4, 
     iter = 1000,
     warmup = 150,
     )

summary(testfit)
plot(testfit, pars = c("mu","base_var","alpha","kappa", "beta1", "beta2"))
plot(testfit, pars = "traits")

samples <- util$extract_expectand_vals(testfit)
dev.off()
par(mfrow=c(3,3))
for(i in 1:9){
  util$plot_expectand_pushforward(samples[[paste0('mu[',i,']')]], 20, "lognormal mean")
  abline(v = mu[i])
}

for(i in 1:9){
  util$plot_expectand_pushforward(samples[[paste0('base_var[',i,']')]], 20, "lognormal var")
  abline(v = base_var[i])
}

par(mfrow = c(2,2))
util$plot_expectand_pushforward(samples[['kappa']], 20, "phylo kappa")
abline(v = kappa)

util$plot_expectand_pushforward(samples[['alpha']], 20, "BA alpha")
abline(v = alpha1)

util$plot_expectand_pushforward(samples[['beta1']], 20, "Beta1 (Competition slope)")
abline(v = Beta1)

util$plot_expectand_pushforward(samples[['beta2']], 20, "Beta2 (on trait distance)")
abline(v = Beta2)
#Heuristic Model
