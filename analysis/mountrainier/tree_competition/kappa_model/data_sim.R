
library(rstan)
library(dplyr)
library(readr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#parameter
set.seed(300)

sigma <- rnorm(1,0, 0.095 / 2.57) #0.05078215
y0 <- rnorm(1,0,10/2.57)#3.354501
k <- rbeta(1,2,2) #0.6317492
r <- rbeta(1,4,2) #0.4923903
beta <- log(1.05)

N <- integer(0)
b <- numeric(0)
bf <- numeric(0)
start_idx <- integer(0)
end_idx <- integer(0)
focal_corr <- numeric(0)

for (i in 1:Nf) {
  id = focal_tags[i]
  this_df <- species_df %>%
    filter(Tag==id)
  
  N <- c(N,length(this_df$"Species(neighbor)"))
  
  #ba <- this_df
  #group_by(`Species(neighbor)`) %>%
  #summarise(basel = sum(ba_n))
  
  b <- c(b,this_df$ba_sum)
  bf<- c(bf,unique(this_df$ba_f))
  
  all_idx <- which(species_df$Tag==id)
  start_idx <- c(start_idx,all_idx[1])
  end_idx <- c(end_idx, all_idx[length(all_idx)])
  
  neighbours <- this_df$`Species(neighbor)`
  focal_spe <- this_df$Species[1]
  for (i in 1:length(neighbours)) {
    neighbour_spe <- neighbours[i]
    if (nchar(focal_spe)==4){
      focal_spe <- toupper(focal_spe)
    }
    if (nchar(neighbour_spe)==4) {
      neighbour_spe <- toupper(neighbour_spe)
    }
    focal_corr <- c(focal_corr,phy_correlation_matrix[focal_spe,neighbour_spe])
  }
}

N_total <- sum(N)

#data input
set.seed(1000)
Nf <- 100
N <- sample(4:8,Nf,replace=TRUE)
N_total <- sum(N)
b <- runif(N_total,2.1,2.5) # x1000cm^2
focal_corr <- runif(N_total,0,1)

# start and end index
start_idx <- numeric(0)
end_idx <- numeric(0)
i <- 1
for (t in N) {
  start_idx <- c(start_idx,i)
  end_idx <- c(end_idx,i+t-1)
  i = i+t
}
bf <- runif(Nf,2.1,2.5)


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
  competition[i]=beta*(baphy[i]-BA_compet0)
  avails[i] = (bf[i]-bf0)*r
  mu[i]=log(y0) + avails[i] - competition[i]
}



log_y <- rnorm(Nf,mu,sigma)
y <- exp(log_y)


stan_data <- list(
  Nf = Nf,
  N = N,
  N_total = N_total,
  b = b,
  focal_corr = focal_corr,
  start_idx = start_idx,
  end_idx = end_idx,
  bf = bf,
  y = y
)

stan_data$Nf <- as.integer(Nf)
stan_data$start_idx <- as.integer(start_idx)
stan_data$end_idx   <- as.integer(end_idx)

#data simulation
fit <- stan(
  file = "model.stan",
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