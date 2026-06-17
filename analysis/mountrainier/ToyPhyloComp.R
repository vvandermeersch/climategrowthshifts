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


#treephylo <- read.tree("analysis/mountrainier/data/rainiertree.tre")
treephylo <- read.tree("data/rainiertree.tre")
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
    ((BAnb[x]/BA0)^-alpha)*Trait.cost(TDnb[x], k=kappa,beta = beta.phylo)}
  )))
}


#Setting up data
tmpneighb <- neighbordata[!is.na(neighbordata$phylodist),]
tmpneighb$BAnb <- (tmpneighb$DBH.neighbor./2)^2
tmpneighb$sp.nums <-  sp.nums[tmpneighb$Species.neighbor.]
tmpneighb <- tmpneighb[,-which(colSums(is.na(tmpneighb))==nrow(tmpneighb))]
tmpneighb <- tmpneighb[-which(rowSums(is.na(tmpneighb))>0),]

#Simulation----
#Simulating phylo effect on basal area
ntrees <- 400
rpt <- 3 #rings per tree

#Simulate base areas of focal species
#Choose random species, sample DBH from observed
spsim <- sample(treephylo$tip.label, size = ntrees, replace = T)
DBHf<- sample(neighbordata.info$DBH[!is.na(neighbordata.info$DBH)], size = ntrees, replace = T)
BAf <- (DBHf/2)^2

Nneighb <- sample(sapply(unique(neighbordata$Treeid), function(x){sum(neighbordata$Treeid==x)}),
                  ntrees, replace = T)

#Simulate latent traits
latenttraits <- fastBM(treephylo)

#Simulate neighbors
species_map <- c(
  TSHE = "Tsuga_heterophylla",
  THPL = "Thuja_plicata",
  ABAM = "Abies_amabilis",
  PSME = "Pseudotsuga_menziesii",
  ABGR = "Abies_grandis",
  TSME = "Tsuga_mertensiana",
  CANO = "Callitropsis_nootkatensis",
  PIMO = "Pinus_monticola",
  ABLA = "Abies_lasiocarpa",
  ACCI = "Abies_lasiocarpa"
)

Neighbors <- list()
valid_tmpneighb <- tmpneighb[!is.na(tmpneighb$Species.neighbor.), ]
for(i in 1:ntrees){
  Neighbors[[i]] <- valid_tmpneighb[sample(1:nrow(valid_tmpneighb), Nneighb[i]),]
  Neighbors[[i]]$Species.neighbor. <- species_map[Neighbors[[i]]$Species.neighbor.]
  for(j in 1:nrow(Neighbors[[i]])){
    sp1 <- spsim[i]
    sp2 <- Neighbors[[i]]$Species.neighbor.[j]
    
    if(!(sp1 %in% names(latenttraits)) || !(sp2 %in% names(latenttraits))){
      stop(paste("Specieanys mismatch:", sp1, sp2))
    }
    
    Neighbors[[i]]$phylodist[j] <- latenttraits[sp1]-latenttraits[sp2]
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
treesp <- rep(as.integer(sp.nums[spsim]), each = rpt)
BaseGrowth <- c(sapply(1:ntrees, function(x){rlnorm(rpt, 
                                                mu[treesp[x]],
                                                base_var[treesp[x]])}))
AllComps <- rep(unlist(sapply(1:ntrees, function(x){
  Comp.fun(Neighbors[[x]]$BAnb, Neighbors[[x]]$phylodist,
           alpha = alpha1, beta.phylo = Beta2, kappa = kappa,
           BA0 = BA0) #*
})), each = rpt)
RingGrowth <- Beta1*BaseGrowth*AllComps #**
hist(RingGrowth, breaks = seq(0,ceiling(max(RingGrowth,na.rm=TRUE)), 0.25))

#Demo rings

#demorings <- read.csv("analysis/mountrainier/data/treerings_ailene/SouthSidecoresX.csv")
demorings <- read.csv("data/treerings_ailene/SouthSidecoresX.csv")
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

start_ids <- rep(c(1, 1+cumsum(Nneighb[-length(Nneighb)])), each = rpt)
end_ids <- rep(cumsum(Nneighb), each = rpt)

standata <- list(
  Nobs = ntrees*rpt, #Number of observations
  S = max(sp.nums), #Number of unique species
  y = RingGrowth, #Ring growth observations
  N_neighbors = sum(Nneighb), #Total number of neighbors
  tree_sp = treesp, #Species of focal tree in each observation
  tree_N_neighbs = rep(Nneighb, each = rpt), #Number of neighbors for each obs
  tree_start_idxs = start_ids, #Neighbor ids for each obs
  tree_end_idxs = end_ids,
  
  neighbor_BA = unlist(c(sapply(1:ntrees, function(x){Neighbors[[x]]$BAnb}))),
  neighbor_sp = unlist(c(sapply(1:ntrees, function(x){sp.nums[Neighbors[[x]]$Species.neighbor.]}))),
  BA0 = 10,
  Sigma = vcv(treephylo)
  )

options(mc.cores = parallel::detectCores())
testfit <- stan("stan/avery_phylo.stan",
     model_name = "testrun!",
     data = standata, 
     chains = 4, 
     iter = 2000,
     warmup = 1000,
     )

summary(testfit)
plot(testfit, pars = c("mu","base_var","alpha","kappa", "beta1", "beta2"))
plot(testfit, pars = "traits")

samples <- util$extract_expectand_vals(testfit)
dev.off()

pdf("posteriro_truth.pdf",width=10,height=8)
on.exit(dev.off())
par(mfrow=c(3,3))
for(i in 1:9){
  util$plot_expectand_pushforward(samples[[paste0('mu[',i,']')]], 20, bquote(mu[.(i)] ~ "(lognormal mean growth)"))
  abline(v = mu[i])
}

for(i in 1:9){
  util$plot_expectand_pushforward(samples[[paste0('base_var[',i,']')]], 20, paste0("mean growth variance[sp",i,"]"))
  abline(v = base_var[i])
}

for(i in 1:9){
  util$plot_expectand_pushforward(samples[[paste0('traits[',i,']')]], 20, paste0("latent traits[sp",i,"]"))
  abline(v = latenttraits[i])
}

par(mfrow = c(2, 2))
util$plot_expectand_pushforward(samples[[paste0('beta1')]], 20, bquote(beta[.(1)]~"(overall competition coefficient)"))
abline(v = Beta1)

util$plot_expectand_pushforward(samples[[paste0('beta2')]], 20, bquote(beta[.(2)]~"(trait-distance decay parameter)"))
abline(v = Beta2)

util$plot_expectand_pushforward(samples[[paste0('alpha')]], 20, bquote(alpha~"(effect of neighbor basal area)"))
abline(v = alpha1)

util$plot_expectand_pushforward(samples[[paste0('kappa')]], 20, bquote(kappa~"(maximum phylogenetic competition strength)"))
abline(v = kappa)

S <- max(sp.nums)
pairwise.dists <- matrix(nrow = S, ncol = S)
pairwise.list <- list()
for(i in 1:(S-1)){
  for(j in (i+1):S){
    pairwise.list[[paste0("pd",i,j)]] <- abs(
      samples[[paste0('traits[',i,']')]] - samples[[paste0('traits[',j,']')]]) 
  }
}

dev.off()
par(mfrow=c(3,3))
for(i in 1:(S-1)){
  #par(mfrow = c(ceiling(sqrt(S-i)), ceiling(sqrt(S-i))))
  for(j in (i+1):S){
  util$plot_expectand_pushforward(pairwise.list[[paste0('pd',i,j)]], 20, paste("distance", i, j))
abline(v = abs(latenttraits[i]-latenttraits[j]))
}}
#Heuristic Model

# prior and posterior and true value
par(mfrow = c(2, 2))

# beta1
prior_beta1 <- rnorm(4000, 1, 1)
plot(density(prior_beta1),
     lty = 2,
     main = "beta1",
     xlim = c(-4, 6),
     ylim = c(0, 0.6))
lines(density(samples[['beta1']]),col='red')
abline(v = Beta1)

# beta2
prior_beta2 <- rnorm(4000, 1, 1)
plot(density(prior_beta2),
     lty = 2,
     main = "beta2")
lines(density(samples[['beta2']]),col='red')
abline(v = Beta2)

# alpha
prior_alpha <- rnorm(4000, 0, 0.1)
plot(density(prior_alpha),
     lty = 2,
     main = "alpha",
     xlim = c(-0.5, 0.5),
     ylim = c(0, 5))
lines(density(samples[['alpha']]),col='red')
abline(v = alpha1)

# kappa
prior_kappa <- rnorm(4000, 0.1, 0.1)
plot(density(prior_kappa),
     lty = 2,
     main = "kappa")
lines(density(samples[['kappa']]),col='red')
abline(v = kappa)

#mu
par(mfrow=c(3,3))

for (i in 1:9){
  prior_mu <- rnorm(4000,0.7,0.3)
  
  d_prior <- density(prior_mu)
  d_post <- density(samples[[paste0('mu[',i,']')]])
  
  plot(
    d_prior,
    lty = 2,
    main = paste0("lognormal mean[", i, "]"),
    xlim = range(c(d_prior$x,d_post$x)),
    ylim = c(0, max(c(d_prior$y, d_post$y)))
  )
  
  lines(d_post,col='red')
  
  abline(v = mu[i])
}

# base var
par(mfrow = c(3,3))

for(i in 1:9){
  
  prior_basevar <- rnorm(4000, 0.5, 0.3)
  
  d_prior <- density(prior_basevar)
  d_post <- density(samples[[paste0('base_var[',i,']')]])
  
  plot(
    d_prior,
    lty = 2,
    main = paste0("base_var[", i, "]"),
    xlim = range(c(d_prior$x,d_post$x)),
    ylim = c(0, max(c(d_prior$y, d_post$y)))
  )
  
  lines(d_post,col='red')
  
  abline(v = base_var[i])
}

# traits
library(MASS)

prior_traits <- MASS::mvrnorm(
  4000,
  mu = rep(0, S),
  Sigma = standata$Sigma
)

par(mfrow = c(3,3))

for(i in 1:9){
  
  d_prior <- density(prior_traits[,i])
  
  d_post <- density(
    samples[[paste0('traits[', i, ']')]]
  )
  
  plot(
    d_prior,
    lty = 2,
    lwd = 2,
    main = paste0("traits[", i, "]"),
    xlim = range(c(d_prior$x,d_post$x)),
    ylim = c(0, max(c(d_prior$y, d_post$y)))
  )
  
  lines(
    d_post,
    col = "red",
    lwd = 2
  )
  
  abline(
    v = latenttraits[i],
    col = "blue",
    lwd = 2
  )
}

# diagnostics
diagnose <- extract_hmc_diagnostics(testfit)
check_all_hmc_diagnostics(diagnose)

expectant_vals <- extract_expectand_vals(testfit)
check_all_expectand_diagnostics(expectant_vals)