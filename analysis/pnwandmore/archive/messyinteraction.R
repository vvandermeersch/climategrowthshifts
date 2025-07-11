
wd <- "~/projects/climategrowthshifts/analysis/pnw"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

fit <- readRDS(file.path(wd, 'output/model', 'fit_28june2025_completepooling_winterac.rds'))
samples <- util$extract_expectand_vals(fit)

alpha <- median(samples[['alpha']])
beta_sm <- median(samples[['beta_sm']])
beta_vpd <- median(samples[['beta_vpd']])
beta_smvpd <- median(samples[['beta_smvpd']])

sm0 <- 25
vpd0 <- 8

sm <- 5
vpd <- 5:30

logrw <- alpha + beta_sm*(sm-sm0) + beta_vpd*(vpd-vpd0) + beta_smvpd *(vpd-vpd0)*(sm-sm0)
rw <- exp(logrw)

plot(rw ~ vpd, col = 'pink', ylim = c(0,2), type = 'l')

sm <- 35
logrw <- alpha + beta_sm*(sm-sm0) + beta_vpd*(vpd-vpd0) + beta_smvpd *(vpd-vpd0)*(sm-sm0)
rw <- exp(logrw)
lines(rw ~ vpd, col = 'purple')


fit <- readRDS(file.path(wd, 'output/model', 'fit_28june2025_completepooling.rds'))
samples <- util$extract_expectand_vals(fit)

alpha <- median(samples[['alpha']])
beta_sm <- median(samples[['beta_sm']])
beta_vpd <- median(samples[['beta_vpd']])

sm0 <- 25
vpd0 <- 8

sm <- 5
vpd <- 5:30

logrw <- alpha + beta_sm*(sm-sm0) + beta_vpd*(vpd-vpd0) 
rw <- exp(logrw)

plot(rw ~ vpd, col = 'pink', ylim = c(0,2), type = 'l')

sm <- 35
logrw <- alpha + beta_sm*(sm-sm0) + beta_vpd*(vpd-vpd0)
rw <- exp(logrw)
lines(rw ~ vpd, col = 'purple')

