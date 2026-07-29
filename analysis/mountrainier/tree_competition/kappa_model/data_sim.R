
library(rstan)
library(dplyr)
library(readr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#parameter
set.seed(300)
S_nu <- 16
sigma <- rnorm(1,0, 0.095 / 2.57) #0.05078215

y0_mu <- rnorm(1,0.9,0.4)
y0_tau <- abs(rnorm(1,0,0.4))
y0_raw <- rnorm(S_nu)

beta_mu <- rnorm(1,-3.8,log(1.8))
beta_tau <- rnorm(1,0,log(1.3))
beta_raw <- rnorm(S_nu)

logodds_mu <- rnorm(1,0.67,0.9)
logodds_tau <- abs(rnorm(1,0,0.1))
logodds_raw <- rnorm(S_nu)

# r_mu <- rbeta(1,8,4)
# r_phi <- rnorm(1,13,2)
# r <- rbeta(S_nu,r_mu*r_phi,(1-r_mu)*r_phi)

k <- rbeta(1,2,2) #0.6317492
#r <- rbeta(S_nu,4,2) #0.4923903
#beta <- log(1.05)
#beta <- rnorm(S_nu,0,log(1.1)/2.57)
#beta <- abs(rnorm(S_nu,0,log(1.1)/2.57));

# data input (real covariates, 6 species)
# N <- integer(0)
# b <- numeric(0)
# ra <- numeric(0)
# start_idx <- integer(0)
# end_idx <- integer(0)
# focal_corr <- numeric(0)
# 
# species_df <- read.csv("data/processed data/species_2008.csv")
# Nf <- length(unique(species_df$Tag))
# for (i in 1:Nf) {
#   id = focal_tags[i]
#   this_df <- species_df %>%
#     filter(Tag==id)
# 
#   N <- c(N,length(this_df$Species.neighbor.))
# 
#   #ba <- this_df
#   #group_by(`Species(neighbor)`) %>%
#   #summarise(basel = sum(ba_n))
# 
#   b <- c(b,this_df$ba_sum)
#   ra<- c(ra,unique(this_df$radius))
# 
#   all_idx <- which(species_df$Tag==id)
#   start_idx <- c(start_idx,all_idx[1])
#   end_idx <- c(end_idx, all_idx[length(all_idx)])
# 
#   neighbours <- this_df$Species.neighbor.
#   focal_spe <- this_df$Species[1]
#   for (i in 1:length(neighbours)) {
#     neighbour_spe <- neighbours[i]
#     if (nchar(focal_spe)==4){
#       focal_spe <- toupper(focal_spe)
#     }
#     if (nchar(neighbour_spe)==4) {
#       neighbour_spe <- toupper(neighbour_spe)
#     }
#     focal_corr <- c(focal_corr,phy_correlation_matrix[focal_spe,neighbour_spe])
#   }
# }
# 
# N_total <- sum(N)
# 
# tree_sp <- species_df %>%
#   select(Tag,Species)%>%
#   unique()%>%
#   pull(Species)
# 
# tree_sp_numeric <- as.numeric(factor(tree_sp))

#data input (fake)
set.seed(1000)
Nf <- 266
N <- sample(4:8,Nf,replace=TRUE)
N_total <- sum(N)
b <- runif(N_total,2.1,2.5) # x1000cm^2
focal_corr <- runif(N_total,0,1)
tree_sp_numeric <- sample(1:S_nu,Nf,replace=TRUE)

# start and end index
start_idx <- numeric(0)
end_idx <- numeric(0)
i <- 1
for (t in N) {
  start_idx <- c(start_idx,i)
  end_idx <- c(end_idx,i+t-1)
  i = i+t
}
ra <- runif(Nf,2,3)


#transformed parameter
BA_compet0 <- 16
#bf0 <- 2
r0 <- 2.5
y0 <- numeric(S_nu)
beta <- numeric(S_nu)
r <- numeric(S_nu)
competition <- numeric(Nf)
baphy <- numeric(Nf)
avails <- numeric(Nf)
mu <- numeric(Nf)
for (i in 1:S_nu) {
  y0[i] <- exp(y0_mu+y0_raw[i]*y0_tau)
  beta[i] <- exp(beta_mu + beta_raw[i]*beta_tau)
  r[i] <- plogis(logodds_mu + logodds_tau*logodds_raw[i])
}

for (i in 1:Nf) {
  bn <- b[start_idx[i]:end_idx[i]]
  corrn <- focal_corr[start_idx[i]:end_idx[i]]
  baphy[i] = sum(bn*(corrn^k))
  competition[i]=beta[tree_sp_numeric[i]]*(baphy[i]-BA_compet0)
  avails[i] = (log(ra[i])-log(r0))*r[tree_sp_numeric[i]]
  mu[i]=log(y0[tree_sp_numeric[i]]) + avails[i] - competition[i]
}


log_y <- rnorm(Nf,mu,sigma)
y <- exp(log_y)


stan_data <- list(
  Nf = Nf,
  S = S_nu,
  N = N,
  tree_sp = tree_sp_numeric,
  N_total = N_total,
  b = b,
  focal_corr = focal_corr,
  start_idx = start_idx,
  end_idx = end_idx,
  ra = ra,
  y = y
)

stan_data$Nf <- as.integer(Nf)
stan_data$start_idx <- as.integer(start_idx)
stan_data$end_idx   <- as.integer(end_idx)

#data simulation
fit <- stan(
  file = 'deltamodel_multispecies.stan',
  data = stan_data,
  iter = 2000,
  chains = 4,
  seed = 123,
)
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

#diagnostics
diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)
samples <- util$extract_expectand_vals(fit)
util$check_all_expectand_diagnostics(samples)

# visualization
par(mfrow = c(1,1))
util$plot_expectand_pushforward(samples[["sigma"]], 25,
                                display_name=bquote(sigma~"(standard deviation of tree growth)"),
                                #baseline=sigma,
                                baseline_col=util$c_mid_teal) 
abline(v=sigma)

util$plot_expectand_pushforward(samples[["k"]],25,display_name = bquote(kappa),
                                baseline_col = util$c_mid_teal)
abline(v=k)

par(mfrow = c(2,3))
for (i in 1:S_nu){
  util$plot_expectand_pushforward(samples[[paste0("y0[",i,"]")]],25,
                                  display_name = paste0("y0[",i,"]"),
                                  baseline_col=util$c_mid_teal)
  abline(v=y0[i])
}

for(i in 1:S_nu){
  util$plot_expectand_pushforward(samples[[paste0("beta[",i,"]")]],25,
                                  display_name = bquote(beta[.(i)]),
                                  baseline_col = util$c_mid_teal)
  abline(v=beta[i])
}

for(i in 1:S_nu){
  util$plot_expectand_pushforward(samples[[paste0("r[",i,"]")]],25,
                                  display_name = paste0("r[",i,"]"),
                                  baseline_col = util$c_mid_teal)
  abline(v=r[i])
}
par(mfrow=c(1,1))
util$plot_expectand_pushforward(samples[["y0_mu"]],25,display_name = bquote(y0_mu),
                                baseline_col = util$c_mid_teal)
abline(v=y0_mu)

util$plot_expectand_pushforward(samples[["y0_tau"]],25,display_name = bquote(y0_tau),
                                baseline_col = util$c_mid_teal)
abline(v=y0_tau)

util$plot_expectand_pushforward(samples[["beta_mu"]],25,display_name = bquote(beta_mu),
                                baseline_col = util$c_mid_teal)
abline(v=beta_mu)

util$plot_expectand_pushforward(samples[["beta_tau"]],25,display_name = bquote(beta_tau),
                                baseline_col = util$c_mid_teal)
abline(v=beta_tau)

util$plot_expectand_pushforward(samples[["r_mu"]],25,display_name = bquote(r_mu),
                                baseline_col = util$c_mid_teal)
abline(v=r_mu)

util$plot_expectand_pushforward(samples[["r_phi"]],25,display_name = bquote(r_phi),
                                baseline_col = util$c_mid_teal)
abline(v=r_phi)
