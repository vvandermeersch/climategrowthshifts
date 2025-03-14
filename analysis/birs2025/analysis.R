######################################################################
#
# Setup
#
######################################################################

par(family="serif", las=1, bty="l",
    cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

######################################################################
#
# Data Exploration
#
######################################################################

raw_data1 <- read.csv('data/tag_E60.csv')
raw_data2 <- read.csv('data/snowfall_longmire.csv')
raw_data3 <- read.csv('data/av06_climate.csv')
raw_data4 <- read.csv('data/snowfreegsl_longmire.csv')

# Ring widths
par(mfrow=c(1, 2))

util$plot_line_hist(raw_data1$rw_core1..mm., 0, 2, 0.1,
                    xlab="Ring Width (mm)", main="Core 1")

util$plot_line_hist(raw_data1$rw_core2..mm., 0, 2, 0.1,
                    xlab="Ring Width (mm)", main="Core 2")

par(mfrow=c(1, 1))

plot(raw_data1$year, raw_data1$rw_core1..mm., 
     pch=16, cex=1.0, xlab="Year", ylab="Ring Width (mm)")
points(raw_data1$year, raw_data1$rw_core2..mm., 
       pch=16, cex=1.0, col=util$c_mid_teal)

# Snowfall
par(mfrow=c(1, 1))

util$plot_line_hist(raw_data2$snowfall..m., 500, 11500, 500,
                    xlab="Snowfall (mm)")

par(mfrow=c(1, 1))

plot(raw_data2$year, raw_data2$snowfall..m., 
     pch=16, cex=1.0, xlab="Year", ylab="Snowfall (mm)")

# GDDB
par(mfrow=c(1, 1))

util$plot_line_hist(raw_data3$gddb5..degC., 
                    1500, 3000, 250,
                    xlab="GDDB (oC)")

par(mfrow=c(1, 1))

plot(raw_data3$year, raw_data3$gddb5..degC., 
     pch=16, cex=1.0, xlab="Year", ylab="GDDB (oC)")

# SWE
par(mfrow=c(1, 1))

util$plot_line_hist(raw_data3$swe..mm., 
                    0, 2500, 250,
                    xlab="SWE (mm)")

par(mfrow=c(1, 1))

plot(raw_data3$year, raw_data3$swe..mm., 
     pch=16, cex=1.0, xlab="Year", ylab="SWE (mm)")

# GSL
par(mfrow=c(1, 1))

util$plot_line_hist(raw_data4$gsl..days., 
                    110, 230, 10,
                    xlab="GSL (days)")

par(mfrow=c(1, 1))

plot(raw_data4$year, raw_data4$gsl..days., 
     pch=16, cex=1.0, xlab="Year", ylab="GSL (days)")


# Format data
data <- list()
data$year <- 1931:2009
data$log_rw_obs <- sapply(data$year, 
                           function(y) 
                             log(raw_data1$rw_core1..mm.[raw_data1$year == y][1]))
data$snowfall <- 1e-3 * sapply(data$year, 
                               function(y) 
                               raw_data2$snowfall..m.[raw_data2$year == y][1])
data$gddb <- sapply(data$year, 
                    function(y) 
                    raw_data3$gddb5..degC.[raw_data3$year == y][1])
data$swe <- sapply(data$year, 
                   function(y) 
                   raw_data3$swe..mm.[raw_data3$year == y][1])
data$gsl <- sapply(data$year, 
                   function(y) 
                     raw_data4$gsl..days.[raw_data4$year == y][1])
data$N <- length(data$year)

######################################################################
#
# Model 1
#
######################################################################

l <- log(0.2)
u <- log(5)

0.5 * (u + l)
0.5 * (u - l) / 2.32

# Prior Checks
fit <- stan(file='stan_programs/model1_prior.stan',
            data=data, seed=5838299,
            warmup=1000, iter=2024, refresh=0)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 'rho', 
                                         'gamma', 'sigma'))
util$check_all_expectand_diagnostics(base_samples)

par(mfrow=c(2, 1))

rw_names <- sapply(1:data$N,
                  function(n) paste0('rw[', n, ']'))
util$plot_realizations(samples, rw_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 10))

util$plot_conn_pushforward_quantiles(samples, rw_names, data$year,
                                     xlab="Year", ylab="Ring Width (mm)", 
                                     display_ylim=c(0, 10))


# Posterior Quantification
fit <- stan(file='stan_programs/model1.stan',
            data=data, seed=5838299,
            warmup=1000, iter=2024, refresh=0)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 'rho', 
                                         'gamma', 'sigma'))
util$check_all_expectand_diagnostics(base_samples)

# Retrodicive Check
par(mfrow=c(1, 1))

rw_names <- sapply(1:data$N,
                   function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conn_pushforward_quantiles(samples, rw_names, data$year,
                                     xlab="Year", ylab="Log Ring Width Per mm", 
                                     display_ylim=c(-2, 2))
points(data$year, data$log_rw_obs, pch=16, cex=1.0, col="white")
points(data$year, data$log_rw_obs, pch=16, cex=0.8, col="black")

# Posterior inferences
par(mfrow=c(2, 1))

rw_names <- sapply(1:data$N,
                   function(n) paste0('rw[', n, ']'))
util$plot_realizations(samples, rw_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")

util$plot_conn_pushforward_quantiles(samples, rw_names, data$year,
                                     xlab="Year", ylab="Ring Width (mm)", 
                                     display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")


par(mfrow=c(2, 2))

util$plot_expectand_pushforward(samples[['sigma']], 40,
                                flim=c(0, 0.27),
                                display_name="sigma")
xs <- seq(-0.15, +0.15, 0.01)
ys <- dnorm(xs, 0, 0.095 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['alpha']], 25,
                                display_name="alpha")
xs <- seq(-1.5, 1.5, 0.01)
ys <- dnorm(xs, 0, 0.69)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['rho']], 25,
                                display_name="rho")
xs <- seq(0, 60, 0.01)
ys <- dlnorm(xs, 3.55, 0.24)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['gamma']], 25,
                                display_name="gamma",)
xs <- seq(0, 3, 0.01)
ys <- dnorm(xs, 0, log(10) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

######################################################################
#
# Model 2
#
######################################################################

l <- log(3)
u <- log(10)

0.5 * (u + l)
0.5 * (u - l) / 2.32

# Posterior Quantification
fit <- stan(file='stan_programs/model2.stan',
            data=data, seed=5838299,
            warmup=1000, iter=2024, refresh=0)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 
                                         'rho', 'gamma', 
                                         'rho_sh', 'gamma_sh',
                                         'sigma'))
util$check_all_expectand_diagnostics(base_samples)

# Retrodicive Check
par(mfrow=c(1, 1))

rw_names <- sapply(1:data$N,
                   function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conn_pushforward_quantiles(samples, rw_names, data$year,
                                     xlab="Year", ylab="Log Ring Width Per mm", 
                                     display_ylim=c(-2, 2))
points(data$year, data$log_rw_obs, pch=16, cex=1.0, col="white")
points(data$year, data$log_rw_obs, pch=16, cex=0.8, col="black")

par(mfrow=c(2, 1))

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$snowfall,
                                       0.5, 11.5, 1, data$log_rw_obs, 
                                       xlab="Snowfall (m)")

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$snowfall,
                                       0.5, 11.5, 1, data$log_rw_obs, 
                                       residual=TRUE, xlab="Snowfall (m)")

par(mfrow=c(2, 1))

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gddb,
                                       1700, 2800, 100, data$log_rw_obs, 
                                       xlab="GDDB (oC)")

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gddb,
                                       1700, 2800, 100, data$log_rw_obs, 
                                       residual=TRUE, xlab="GDDB (oC)")

par(mfrow=c(2, 1))

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$swe,
                                       0, 2250, 250, data$log_rw_obs, 
                                       xlab="SWE (mm)")

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$swe,
                                       0, 2250, 250, data$log_rw_obs, 
                                       residual=TRUE, xlab="SWE (mm)")

par(mfrow=c(2, 1))

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gsl,
                                       110, 230, 5, data$log_rw_obs, 
                                       xlab="GSL (days)")

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gsl,
                                       110, 230, 5, data$log_rw_obs, 
                                       residual=TRUE, xlab="GSL (days)")



# Posterior inferences
par(mfrow=c(2, 1))

rw_names <- sapply(1:data$N,
                   function(n) paste0('rw[', n, ']'))
util$plot_realizations(samples, rw_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")

util$plot_conn_pushforward_quantiles(samples, rw_names, data$year,
                                     xlab="Year", ylab="Ring Width (mm)", 
                                     display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")


par(mfrow=c(3, 2))

util$plot_expectand_pushforward(samples[['sigma']], 40, 
                                flim=c(0, 0.22),
                                display_name="sigma")
xs <- seq(-0.15, +0.15, 0.01)
ys <- dnorm(xs, 0, 0.095 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['alpha']], 25,
                                display_name="alpha")
xs <- seq(-1.5, 1.5, 0.01)
ys <- dnorm(xs, 0, 0.69)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['rho']], 25,
                                display_name="rho")
xs <- seq(0, 60, 0.01)
ys <- dlnorm(xs, 3.55, 0.24)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['gamma']], 25,
                                display_name="gamma")
xs <- seq(0, 3, 0.01)
ys <- dnorm(xs, 0, log(10) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['rho_sh']], 25,
                                display_name="rho_sh")
xs <- seq(0, 60, 0.01)
ys <- dlnorm(xs, 1.7, 0.26)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['gamma_sh']], 25,
                                display_name="gamma_sh")
xs <- seq(0, 3, 0.01)
ys <- dnorm(xs, 0, log(3) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)


######################################################################
#
# Model 3
#
######################################################################

g <- function(x, xo, kappa) {
  x**kappa / (xo + x**kappa)
}

xo <- 2
kappa <- 100

xs <- seq(0, 5, 0.01)
ys <- sapply(xs, function(x) g(x, xo, kappa))
plot(xs, ys, type="l", lwd=2)


f <- function(x, beta, gamma, kappa) {
  f1 <- (exp(-beta * x) + gamma)
  f2 <- x / (beta * kappa + x)
  f1 * f2
}

beta <- 5
gamma <- 0.04
kappa <- 0.05

par(mfrow=c(1, 1))

xs <- seq(0, 5, 0.01)
ys <- sapply(xs, function(x) f(x, beta, gamma, kappa))
plot(xs, ys, type="l", lwd=2)
abline(h=gamma, lwd=2, lty=2, col="#DDDDDD")

# Posterior Quantification
fit <- stan(file='stan_programs/model3.stan',
            data=data, seed=5838299,
            warmup=1000, iter=2024, refresh=0)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 'beta_gsl', 
                                         'rho', 'gamma', 
                                         'rho_sh', 'gamma_sh',
                                         'sigma'))
util$check_all_expectand_diagnostics(base_samples)

# Retrodicive Check
par(mfrow=c(1, 1))

rw_names <- sapply(1:data$N,
                   function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conn_pushforward_quantiles(samples, rw_names, data$year,
                                     xlab="Year", ylab="Log Ring Width Per mm", 
                                     display_ylim=c(-2, 2))
points(data$year, data$log_rw_obs, pch=16, cex=1.0, col="white")
points(data$year, data$log_rw_obs, pch=16, cex=0.8, col="black")

par(mfrow=c(2, 1))

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$snowfall,
                                       0.5, 11.5, 1, data$log_rw_obs, 
                                       xlab="Snowfall (m)")

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$snowfall,
                                       0.5, 11.5, 1, data$log_rw_obs, 
                                       residual=TRUE, xlab="Snowfall (m)")

par(mfrow=c(2, 1))

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gsl,
                                       110, 230, 5, data$log_rw_obs, 
                                       xlab="GSL (days)")

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gsl,
                                       110, 230, 5, data$log_rw_obs, 
                                       residual=TRUE, xlab="GSL (days)")

# Posterior inferences
par(mfrow=c(3, 1))

mu_names <- sapply(1:data$N,
                   function(n) paste0('mu0[', n, ']'))
util$plot_realizations(samples, mu_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")

mu_names <- sapply(1:data$N,
                   function(n) paste0('mu1[', n, ']'))
util$plot_realizations(samples, mu_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")

rw_names <- sapply(1:data$N,
                   function(n) paste0('rw[', n, ']'))
util$plot_realizations(samples, rw_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")

par(mfrow=c(1, 1))

rw_names <- sapply(1:data$N,
                   function(n) paste0('rw[', n, ']'))
util$plot_conn_pushforward_quantiles(samples, rw_names, data$year,
                                     xlab="Year", ylab="Ring Width (mm)", 
                                     display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")


par(mfrow=c(1, 3))

util$plot_expectand_pushforward(samples[['sigma']], 40, 
                                flim=c(0, 0.22),
                                display_name="sigma")
xs <- seq(-0.15, +0.15, 0.01)
ys <- dnorm(xs, 0, 0.095 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['alpha']], 25,
                                display_name="alpha")
xs <- seq(-1.5, 1.5, 0.01)
ys <- dnorm(xs, 0, 0.69)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_gsl']], 25,
                                display_name="beta_gsl")
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

par(mfrow=c(2, 2))

util$plot_expectand_pushforward(samples[['rho']], 25,
                                display_name="rho")
xs <- seq(0, 60, 0.01)
ys <- dlnorm(xs, 3.55, 0.24)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['gamma']], 25,
                                display_name="gamma")
xs <- seq(0, 3, 0.01)
ys <- dnorm(xs, 0, log(10) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['rho_sh']], 25,
                                display_name="rho_sh")
xs <- seq(0, 60, 0.01)
ys <- dlnorm(xs, 1.7, 0.26)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['gamma_sh']], 25,
                                display_name="gamma_sh")
xs <- seq(0, 3, 0.01)
ys <- dnorm(xs, 0, log(3) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)


######################################################################
#
# Data Exploration
#
######################################################################

raw_data1 <- read.csv('data/abam_AV06.csv')
raw_data2 <- read.csv('data/snowfreegsl_longmire.csv')

# Average ring widths when both cores are present
robust_ave <- function(x, y) {
  ifelse(is.na(y), x, mean(c(x, y)))
}

raw_data1$rw_core_ave..mm. <- as.vector(by(raw_data1, 
                                           seq_len(nrow(raw_data1)), 
                                           function(r) robust_ave(r[[4]], r[[5]])))

uniq_tree_ids <- unique(raw_data1$tree_id)
N_trees <- length(uniq_tree_ids)

# Plot ring width data
par(mfrow=c(5, 2))

for (id in uniq_tree_ids) {
  rw_ave <- raw_data1$rw_core_ave..mm.[raw_data1$tree_id == id]
  util$plot_line_hist(rw_ave, 0, 3.5, 0.1,
                    xlab="Ring Width (mm)", 
                    main=paste("Core Average, Tree", id))
}

par(mfrow=c(5, 2))

for (id in uniq_tree_ids) {
  year <- raw_data1$year[raw_data1$tree_id == id]
  rw_ave <- raw_data1$rw_core_ave..mm.[raw_data1$tree_id == id]
  
  plot(year, rw_ave, pch=16, cex=1.0, 
       xlab="Year", ylab="Ring Width (mm)", ylim=c(0, 3.5),
       main=paste("Tree", id))
}

# Plot GSL data
par(mfrow=c(2, 1))

util$plot_line_hist(raw_data4$gsl..days., 
                    110, 230, 10,
                    xlab="GSL (days)")

plot(raw_data4$year, raw_data4$gsl..days., 
     pch=16, cex=1.0, xlab="Year", ylab="GSL (days)")

# Format data into ragged arrays
all_years <- intersect(raw_data1$year, raw_data2$year)
N_all_years <- length(all_years)

uniq_tree_ids <- unique(raw_data1$tree_id)
N_trees <- length(uniq_tree_ids)

log_rw_obs <- c()
gsl_obs <- c()
years <- c()
all_years_idxs <- c()
N_years <- c()

idx <- 1
tree_start_idxs <- c()
tree_end_idxs <- c()

for (id in uniq_tree_ids) {
  raw_data1_tree <- raw_data1[raw_data1$tree_id == id,]
  
  years_tree <- intersect(raw_data1_tree$year, raw_data2$year)
  all_years_idxs_tree <- sapply(years_tree, function(y) which(all_years == y))
  N_years_tree <- length(years_tree)
  
  log_rw_obs_tree <- sapply(years_tree, 
                            function(y) 
                            log(raw_data1_tree$rw_core_ave..mm.[raw_data1_tree$year == y][1]))
  gsl_obs_tree <- sapply(years_tree, 
                         function(y) 
                         raw_data2$gsl..days.[raw_data2$year == y][1])
  
  log_rw_obs <- c(log_rw_obs, log_rw_obs_tree)
  gsl_obs <- c(gsl_obs, gsl_obs_tree)
  years <- c(years, years_tree)
  all_years_idxs <- c(all_years_idxs, all_years_idxs_tree)
  N_years <- c(N_years, N_years_tree)
  
  tree_start_idxs <- c(tree_start_idxs, idx)
  idx <- idx + N_years_tree
  tree_end_idxs <- c(tree_end_idxs, idx - 1)
}

# Cross check sizes
N_trees

length(log_rw_obs)
length(gsl_obs)
length(years)
length(all_years_idxs)
length(N_years)

length(tree_start_idxs)
length(tree_end_idxs)

# Plot all trees using new data format
par(mfrow=c(5, 2))

for (t in 1:N_trees) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]

  year <- years[idxs]
  rw_ave <- exp(log_rw_obs[idxs])

  plot(year, rw_ave, pch=16, cex=1.0, 
       xlab="Year", ylab="Ring Width Per mm", ylim=c(0, 3.5),
       main=paste("Tree", uniq_tree_ids[t]))
}


# Cross check
par(mfrow=c(2, 1))

t <- 12

id <- uniq_tree_ids[t]
year <- raw_data1$year[raw_data1$tree_id == id]
rw_ave <- raw_data1$rw_core_ave..mm.[raw_data1$tree_id == id]

plot(year[year >= 1931], rw_ave[year >= 1931], pch=16, cex=1.0, 
     xlab="Year", ylab="Ring Width (mm)", ylim=c(0, 3.5),
     main=paste("Tree", uniq_tree_ids[t]))

idxs <- tree_start_idxs[t]:tree_end_idxs[t]
year <- years[idxs]
rw_ave <- exp(log_rw_obs[idxs])

plot(year, rw_ave, pch=16, cex=1.0, 
     xlab="Year", ylab="Ring Width Per mm", ylim=c(0, 3.5),
     main=paste("Tree", uniq_tree_ids[t]))

# Collection data into list
N <- length(years)

data <- mget(c('N', 'N_all_years', 'N_trees', 
               'log_rw_obs', 'gsl_obs', 
               'all_years', 'years', 'all_years_idxs', 'N_years', 
               'tree_start_idxs', 'tree_end_idxs'))

# Posterior Quantification
fit <- stan(file='stan_programs/model4.stan',
            data=data, seed=5838299,
            warmup=1000, iter=2024, refresh=10)

diagnostics <- util$extract_hmc_diagnostics(fit)
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)
base_samples <- util$filter_expectands(samples,
                                       c('alpha', 'beta_gsl', 
                                         'rho', 'gamma',
                                         'f_tilde_sh',
                                         'rho_sh', 'gamma_sh',
                                         'sigma'),
                                         check_arrays=TRUE)
util$check_all_expectand_diagnostics(base_samples)

# Retrodicive Check
par(mfrow=c(5, 2))

for (t in 1:data$N_trees) {
  idxs <- tree_start_idxs[t]:tree_end_idxs[t]
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(samples, rw_names, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-3, 3))
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
}

par(mfrow=c(2, 1))

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gsl,
                                       110, 230, 5, data$log_rw_obs, 
                                       xlab="GSL (days)")

pred_names <- sapply(1:data2$N, function(n) paste0('log_rw_pred[', n, ']'))
util$plot_conditional_median_quantiles(samples, pred_names, data$gsl,
                                       110, 230, 5, data$log_rw_obs, 
                                       residual=TRUE, xlab="GSL (days)")

# Posterior inferences
par(mfrow=c(1, 1))

rw_names <- sapply(1:data$N_all_years,
                   function(n) paste0('f_sh[', n, ']'))
util$plot_realizations(samples, rw_names, data$all_years,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(-3, 3))
abline(v=1948, lwd=2, lty=2, col="#DDDDDD")
abline(v=1966, lwd=2, lty=2, col="#DDDDDD")
abline(v=1974, lwd=2, lty=2, col="#DDDDDD")
abline(v=2001, lwd=2, lty=2, col="#DDDDDD")

par(mfrow=c(3, 1))

mu_names <- sapply(1:data$N,
                   function(n) paste0('mu0[', n, ']'))
util$plot_realizations(samples, mu_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")

mu_names <- sapply(1:data$N,
                   function(n) paste0('mu1[', n, ']'))
util$plot_realizations(samples, mu_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")

rw_names <- sapply(1:data$N,
                   function(n) paste0('rw[', n, ']'))
util$plot_realizations(samples, rw_names, data$year,
                       xlab="Year", ylab="Ring Width (mm)", 
                       display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")

par(mfrow=c(1, 1))

rw_names <- sapply(1:data$N,
                   function(n) paste0('rw[', n, ']'))
util$plot_conn_pushforward_quantiles(samples, rw_names, data$year,
                                     xlab="Year", ylab="Ring Width (mm)", 
                                     display_ylim=c(0, 2))
points(data$year, exp(data$log_rw_obs), pch=16, cex=1.0, col="white")
points(data$year, exp(data$log_rw_obs), pch=16, cex=0.8, col="black")


par(mfrow=c(1, 3))

util$plot_expectand_pushforward(samples[['sigma']], 40, 
                                flim=c(0, 0.22),
                                display_name="sigma")
xs <- seq(-0.15, +0.15, 0.01)
ys <- dnorm(xs, 0, 0.095 / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['alpha']], 25,
                                display_name="alpha")
xs <- seq(-1.5, 1.5, 0.01)
ys <- dnorm(xs, 0, 0.69)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['beta_gsl']], 25,
                                display_name="beta_gsl")
xs <- seq(0, 1, 0.01)
ys <- dnorm(xs, 0, log(1.8) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

par(mfrow=c(2, 2))

util$plot_expectand_pushforward(samples[['rho']], 25,
                                display_name="rho")
xs <- seq(0, 60, 0.01)
ys <- dlnorm(xs, 3.55, 0.24)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['gamma']], 25,
                                display_name="gamma")
xs <- seq(0, 3, 0.01)
ys <- dnorm(xs, 0, log(10) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['rho_sh']], 25,
                                display_name="rho_sh")
xs <- seq(0, 60, 0.01)
ys <- dlnorm(xs, 1.7, 0.26)
lines(xs, ys, lwd=2, col=util$c_light)

util$plot_expectand_pushforward(samples[['gamma_sh']], 25,
                                display_name="gamma_sh")
xs <- seq(0, 3, 0.01)
ys <- dnorm(xs, 0, log(3) / 2.57)
lines(xs, ys, lwd=2, col=util$c_light)
