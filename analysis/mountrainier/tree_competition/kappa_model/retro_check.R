
library(rstan)
library(dplyr)
library(readr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# posterior
env1 <- new.env()
source('confronting_with_real_data.R',local=env1)
postfit <- env1$fit
post <- rstan::extract(postfit)

# data input
Nf <- env1$Nf
S <- env1$S
N <- env1$N
tree_sp <- env1$tree_sp_numeric
N_total <- env1$N_total
b <- env1$b
focal_corr <- env1$focal_corr
start_idx <- env1$start_idx
end_idx <- env1$end_idx
bf <- env1$bf
y_obs <- env1$y

#parameter
set.seed(300)
draw <- sample(1:length(post$k),6)
y_retro <- matrix(NA,
                  nrow = Nf,
                  ncol = length(draw))
for (t in 1:length(draw)) {
  sigma <- post$sigma[draw[t]]
  k <- post$k[draw[t]]
  y0 <- post$y0[draw[t],]
  beta <- post$beta[draw[t],]
  r <- post$r[draw[t],]

  #transformed parameter
  BA_compet0 <- 16
  bf0 <- 2
  competition <- numeric(Nf)
  baphy <- numeric(Nf)
  avails <- numeric(Nf)
  mu <- numeric(Nf)
  
  
  for (i in 1:Nf) {
    bn <- b[start_idx[i]:end_idx[i]]
    corrn <- focal_corr[start_idx[i]:end_idx[i]]
    baphy[i] = sum(bn*(corrn^k))
    competition[i]=beta[tree_sp[i]]*(baphy[i]-BA_compet0)
    avails[i] = (bf[i]-bf0)*r[tree_sp[i]]
    mu[i]=log(y0[tree_sp[i]]) + avails[i] - competition[i]
  }
  
  log_y <- rnorm(Nf,mu,sigma)
  y_retro[,t] <- exp(log_y)
}
xlim <- range(c(y_obs, y_retro), na.rm = TRUE)
breaks <- seq(xlim[1], xlim[2], length.out = 20)

source("mcmc_visualization_tools.R")
util$plot_hist_quantiles()
# par(mfrow = c(1,2))
# for (c in 1:length(draw)) {
#   hist(y_retro[, c],
#        probability = TRUE,
#        col = "grey",
#        border = "white",
#        xlim = xlim,
#        breaks = breaks,
#        main=paste0("Posterior Draw",c))
#   
#   hist(y_obs,
#        probability = TRUE,
#        add = TRUE,
#        breaks = breaks,
#        col = adjustcolor("red",alpha.f=0.2),
#        border = adjustcolor("red",alpha.f=0.2),
#        lwd = 2)
# }
