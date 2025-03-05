
library(ggplot2)

# We got our model
gompertzmodel <- stan_model(file.path(wd, "stan/heterogeneous_gompertz.stan"))

# Let's simulate data (post-model, pre-data)
gompertz <- function(t, d_max, beta, delta_t_half) {
    return(d_max * exp(-log(2) * exp( -   (2 / log(2)) * (beta / d_max) * (t - (2000 + delta_t_half)))))
}

N_trees <- 20
Nobs_pertree <- round(runif(N_trees,4,10))
Sobs <- c(0,cumsum(Nobs_pertree))
ID_trees <- rep(1:N_trees, times = Nobs_pertree)

years_init <- sample(1970:1980, N_trees, replace = TRUE)
years <- c()
for(t in 1:N_trees){
  yt <- c(years_init[t])
  for(n in 2:Nobs_pertree[t]){
    ytn <- yt[n-1] + sample(4:8, 1)
    yt <- c(yt, ytn )
  }
  years <- c(years, yt)
}

d0 <- runif(N_trees, 10, 60) 
delta_d <- runif(N_trees, 10, 60) # Difference between minimum and maximum DBH (cm)
beta <- runif(N_trees, 0.1, 5) # Intermediate linear growth rate (cm / year)
delta_t_half <- runif(N_trees, -15, 15)

yhat <- c()
for(i in 1:max(Sobs)){
  yhati <- gompertz(years[i], delta_d[ID_trees[i]], beta[ID_trees[i]], delta_t_half[ID_trees[i]]) + d0[ID_trees[i]]
  yhat <- c(yhat, yhati)
}
y <- rnorm(max(Sobs), yhat, 2);

ggplot() +
  geom_point(aes(x = years, y = y)) +
  geom_line(aes(x = years, y = yhat, group = ID_trees)) +
  theme_bw()


# Now, fit simulated data
mdl.data <- list(N = max(Sobs),
                 N_trees =  N_trees,
                 Sobs = Sobs,
                 years = years,
                 dbhs = y)

fit <- sampling(gompertzmodel, mdl.data, 
                chains = 1, 
                iter = 4000, warmup = 3000,
                verbose = TRUE)
