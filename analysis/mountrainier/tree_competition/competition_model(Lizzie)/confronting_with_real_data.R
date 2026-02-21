library(rstan)
library(dplyr)
library(readr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

species_df <- read_csv("data/processed data/species_2008.csv")
phy_correlation_matrix <- readRDS("phy_correlation_matrix.rds")

# data input
set.seed(1000)
focal_tags <- unique(species_df$Tag)
Nf <- length(focal_tags)
#N <- sample(4:8,Nf,replace=TRUE)

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
  
  N <- c(N,length(this_df$`Species(neighbor)`))
  
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


# # transformed parameter
# BA_compet0 <- 16
# bf0 <- 2
# competition <- numeric(Nf)
# baphy <- numeric(Nf)
# avails <- numeric(Nf)
# mu <- numeric(Nf)
# 
# for (i in 1:Nf) {
#   bn <- b[start_idx[i]:end_idx[i]]
#   corrn <- focal_corr[start_idx[i]:end_idx[i]]
#   baphy[i] = sum(bn*(corrn^k))
#   competition[i]=beta*(baphy[i]-BA_compet0)
#   avails[i] = (bf[i]-bf0)*r
#   mu[i]=log(y0) + avails[i] - competition[i]
# }
# 
# 
 focal_unique <- unique(species_df[,c("Tag","2008")])
 y <- focal_unique$`2008`

# stan data
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

fit <- stan(
  file = "model.stan",
  data = stan_data,
  iter = 2000,
  chains = 4,
  seed = 123,
)