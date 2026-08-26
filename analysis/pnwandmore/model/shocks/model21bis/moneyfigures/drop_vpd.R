rm(list = ls())
library(cmdstanr)
library(doFuture)
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('mcmc_custom_functions.R', local = util)
setwd(wd)


fit <- readRDS(file.path(wd, 'output/model/model21bis', 'fit_model21bis_HGSP_24species_365stands_withinit_4threads_climateshocks.rds'))
data <- readRDS(file.path(wd, 'output/model' ,'data_10july2026_24species_365stands_19502022.rds'))

params <- c(
  "mu_alpha", "sigma_alpha", "alpha",
  "sigma_alpha_stand", "alpha_stand",
  "mu_beta_gdd", "sigma_beta_gdd", "beta_gdd",
  "mu_beta_pre", "sigma_beta_pre", "beta_pre",
  "mu_beta_vpd", "sigma_beta_vpd", "beta_vpd",
  "tau_clim", 
  'log_kappa_clim', 'kappa_clim',
  "mu_log_tau_conc", "sigma_log_tau_conc", "tau_conc",
  "mu_phi", "sigma_phi", "logit_phi_sck", "beta_phi_vpd", "beta_phi_pre",
  "mu_omega_conc", "sigma_omega_conc", "logit_omega_conc_sck",
  "mu_omega_shutdown", "sigma_omega_shutdown", "logit_omega_shutdown",
  "sigma"
)

base_samples <- fit$draws(params)
names <- dimnames(base_samples)$variable
base_samples <- lapply(1:dim(base_samples)[3],
                       function(k){t(matrix(base_samples[1:dim(base_samples)[1],1:dim(base_samples)[2],k],
                                            nrow = dim(base_samples)[1], ncol = dim(base_samples)[2]))})
names(base_samples) <- names

probs <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
gdd0 = 10
pre0 = 5
vpd0 = 23
g <- mean(data$gdd_obs)
p <- mean(data$pre_obs)
quantile(data$vpd_obs, c(0.01, 0.99))
vpds <- seq(8, 36, 1)
ntrees <- 1e4

plan(multisession, workers = 20)
registerDoFuture()

output <- foreach(v = vpds, .options.future = list(seed = TRUE)) %dofuture% {
  
  mu <- base_samples[['mu_alpha']] +
    base_samples[['mu_beta_gdd']] * (g - gdd0) +
    base_samples[['mu_beta_pre']] * (p - pre0) +
    base_samples[['mu_beta_vpd']] * (v - vpd0)
  
  phi_sck <- boot::inv.logit(
    base_samples[['mu_phi']] +
      base_samples[['beta_phi_vpd']] * (v - vpd0) +
      base_samples[['beta_phi_pre']] * (p - pre0)
  )
  
  omega <- boot::inv.logit(base_samples[['mu_omega_conc']])
  omega_shutdown <- boot::inv.logit(base_samples[['mu_omega_shutdown']])
  
  tau_conc <- exp(base_samples[['mu_log_tau_conc']])
  
  rw_mat <- matrix(NA, nrow = nrow(mu), ncol = ncol(mu) * ntrees)
  
  for (r in 1:nrow(mu)){
    for (c in 1:ncol(mu)){
      
      for (n in 1:ntrees){
        stand_event <- rbinom(1, 1, phi_sck[r, c])
        if (stand_event){
          tree_shock <- rbinom(1, 1, omega[r, c])
          if (tree_shock){
            shutdown <- rbinom(1, 1, omega_shutdown[r, c])
            if (shutdown){
              rw <- 0
            }else{
              delta_shock <- abs(rnorm(1, 0, tau_conc[r, c]))
              rw <- exp(mu[r, c] - delta_shock)
            }
          }else{
            rw <- exp(mu[r, c])
          }
        }else{
          rw <- exp(mu[r, c])
        }
        rw_mat[r, (c - 1) * ntrees + n] <- rw
      }
      # rw_mat[r, c] <- mean(rw_trees)
    }
  }
  
  list(
    phi = util$ensemble_mcmc_quantile_est(phi_sck, probs),
    mu = util$ensemble_mcmc_quantile_est(mu, probs),
    expmu = util$ensemble_mcmc_quantile_est(exp(mu), probs),
    logrw = util$ensemble_mcmc_quantile_est(log(rw_mat), probs),
    rw = util$ensemble_mcmc_quantile_est(rw_mat, probs)
  )
}
plan(sequential)

phis <- sapply(output, function(l){l[["phi"]]})
rws <- sapply(output, function(l){l[["rw"]]})
logrws <- sapply(output, function(l){l[["logrw"]]})
mus <- sapply(output, function(l){l[["mu"]]})
expmus <- sapply(output, function(l){l[["expmu"]]})

# par(mfrow = c(1,2), mar = c(4,4,1,1))
# plot(y = NULL, x = NULL, xlab = 'VPD (hPa)', ylab = 'ring width (mm)', bty = 'n',
#      xlim = c(10, 50), ylim =  c(0,1.5))
# 
# polygon(c(vpds, rev(vpds)), c(expmus['5%', ], rev(expmus['95%',])),
#         col = '#dae2e2', border = NA)
# polygon(c(vpds, rev(vpds)), c(rws['5%', ], rev(rws['95%',])),
#         col = "#eedede", border = NA)
# 
# 
# lines(x = vpds, y = rws['95%',], col = util$c_mid, lty = 2)
# lines(x = vpds, y = rws['5%',], col = util$c_mid, lty = 2)
# lines(x = vpds, y = rws['50%',], col = util$c_mid)
# 
# lines(x = vpds, y = expmus['95%',], col = util$c_mid_teal, lty = 2)
# lines(x = vpds, y = expmus['5%',], col = util$c_mid_teal, lty = 2)
# lines(x = vpds, y = expmus['50%',], col = util$c_mid_teal)
# 
# plot(y = NULL, x = NULL, xlab = 'VPD (hPa)', ylab = 'p(stand extreme event)', bty = 'n',
#      xlim = c(1, 50), ylim =  c(0,1))
# 
# polygon(c(vpds, rev(vpds)), c(phis['5%',], rev(phis['95%',])),
#         col = "#eedede", border = NA)
# lines(x = vpds, y = phis['95%',], col = util$c_mid, lty = 2)
# lines(x = vpds, y = phis['5%',], col = util$c_mid, lty = 2)
# lines(x = vpds, y = phis['50%',], col = util$c_mid)
# 
# 
# par(mfrow = c(1,2), mar = c(4,4,1,1))
# 
# plot(y = NULL, x = NULL, xlab = 'VPD (hPa)', ylab = 'ring width (mm)', bty = 'n',
#      xlim = c(10, 50), ylim =  c(0,1.25))
# 
# polygon(c(vpds, rev(vpds)), c(expmus['5%', ], rev(expmus['95%',])),
#         col = '#dae2e2', border = NA)
# polygon(c(vpds, rev(vpds)), c(rws['5%', ], rev(rws['95%',])),
#         col = "#eedede", border = NA)
# 
# lines(x = vpds, y = expmus['95%',], col = util$c_mid_teal, lty = 2)
# lines(x = vpds, y = expmus['5%',], col = util$c_mid_teal, lty = 2)
# lines(x = vpds, y = expmus['50%',], col = util$c_mid_teal)
# 
# lines(x = vpds, y = rws['95%',], col = util$c_mid, lty = 2)
# lines(x = vpds, y = rws['5%',], col = util$c_mid, lty = 2)
# lines(x = vpds, y = rws['50%',], col = util$c_mid)
# 
# 
# plot(y = NULL, x = NULL, xlab = 'VPD (hPa)', ylab = 'log(ring width)', bty = 'n',
#      xlim = c(10, 50), ylim =  c(-2,0.5))
# 
# polygon(c(vpds, rev(vpds)), c(mus['5%', ], rev(mus['95%',])),
#         col = '#dae2e2', border = NA)
# polygon(c(vpds, rev(vpds)), c(logrws['5%', ], rev(logrws['95%',])),
#         col = "#eedede", border = NA)
# 
# lines(x = vpds, y = mus['95%',], col = util$c_mid_teal, lty = 2)
# lines(x = vpds, y = mus['5%',], col = util$c_mid_teal, lty = 2)
# lines(x = vpds, y = mus['50%',], col = util$c_mid_teal)
# 
# lines(x = vpds, y = logrws['95%',], col = util$c_mid, lty = 2)
# lines(x = vpds, y = logrws['5%',], col = util$c_mid, lty = 2)
# lines(x = vpds, y = logrws['50%',], col = util$c_mid)
# 
# par(mfrow = c(1,1), mar = c(4,4,1,1))
# 
# plot(y = NULL, x = NULL, xlab = 'VPD (hPa)', ylab = 'ring width (mm)', bty = 'n',
#      xlim = c(0, 50), ylim =  c(0,1.5))
# 
# # polygon(c(vpds, rev(vpds)), c(expmus['5%', ], rev(expmus['95%',])),
# #         col = '#dae2e2', border = NA)
# polygon(c(vpds, rev(vpds)), c(rws['5%', ], rev(expmus['5%',])),
#         col = "#eedede", border = NA)
# polygon(c(vpds, rev(vpds)), c(rws['5%', ]-0.1, rev(expmus['5%',]+0.1)),
#         col = "white", border = NA, 
#         density = 5, angle = 45, lwd = 10)
# 
# lines(x = vpds, y = expmus['95%',], col = util$c_mid_teal, lty = 2, lwd = 1.5)
# lines(x = vpds, y = expmus['5%',], col = util$c_mid_teal, lty = 2, lwd = 1.5)
# lines(x = vpds, y = expmus['50%',], col = util$c_mid_teal, lwd = 1.5)
# 
# lines(x = vpds, y = rws['95%',], col = util$c_mid, lty = 2, lwd = 1.5)
# lines(x = vpds, y = rws['5%',], col = util$c_mid, lty = 2, lwd = 1.5)
# lines(x = vpds, y = rws['50%',], col = util$c_mid, lwd = 1.5)
# 
# 
# plot(y = NULL, x = NULL, xlab = 'VPD (hPa)', ylab = 'ring width (mm)', bty = 'n',
#      xlim = c(1, 50), ylim =  c(0,1.5))
# 
# 
# polygon(c(vpds, rev(vpds)), c(expmus['5%', ], rev(rws['95%',])),
#         col = "#f4ded7", border = NA)
# 
# polygon(c(vpds, rev(vpds)), c(expmus['5%', ], rev(rws['95%',])),
#         col = "#d5e5e2", border = NA, 
#         density = 4, angle = 45, lwd = 10)
# 
# polygon(c(vpds, rev(vpds)), c(rws['95%', ], rev(rws['95%',]+0.1)),
#         col = "white", border = NA)
# polygon(c(vpds, rev(vpds)), c(rws['5%', ], rev(rws['5%',]-0.1)),
#         col = "white", border = NA)
# rect(xleft = 0, xright = 1, ybottom = 0.5, ytop = 2, col = 'white', border = NA)
# rect(xleft = 50, xright = 51, ybottom = 0, ytop = 2, col = 'white', border = NA)
# 
# polygon(c(vpds, rev(vpds)), c(rws['95%', ], rev(expmus['95%',])),
#         col = "#d5e5e2", border = NA)
# polygon(c(vpds, rev(vpds)), c(expmus['5%', ], rev(rws['5%',])),
#         col = "#f4ded7", border = NA)
# 
# lines(x = vpds, y = expmus['95%',], col = "#2F7D6D", lty = 2, lwd = 1.5)
# lines(x = vpds, y = expmus['5%',], col = "#2F7D6D", lty = 2, lwd = 1.5)
# lines(x = vpds, y = expmus['50%',], col = "#2F7D6D", lwd = 2)
# 
# lines(x = vpds, y = rws['95%',], col = "#C75B39", lty = 2, lwd = 1.5)
# lines(x = vpds, y = rws['5%',], col = "#C75B39", lty = 2, lwd = 1.5)
# lines(x = vpds, y = rws['50%',], col = "#C75B39", lwd = 2)


par(mfrow=c(1,1))
plot(y = NULL, x = NULL,
     xlab = "Summer VPD (hPa)", ylab = "Ring width (mm)",
     bty = "n",
     xlim = c(8, 36), ylim = c(0, 1.5))

polygon(c(vpds, rev(vpds)),
        c(expmus["5%", ], rev(expmus["95%", ])),
        col = "#2F7D6D33", border = NA)

polygon(c(vpds, rev(vpds)),
        c(rws["5%", ], rev(rws["95%", ])),
        col = "#C75B3933", border = NA)


# Expected growth
lines(vpds, expmus["95%", ], col = "#2F7D6D", lty = 2, lwd = 1.2)
lines(vpds, expmus["5%",  ], col = "#2F7D6D", lty = 2, lwd = 1.2)
lines(vpds, expmus["50%", ], col = "#2F7D6D", lwd = 2.5)

# Ring width
lines(vpds, rws["95%", ], col = "#C75B39", lty = 2, lwd = 1.2)
lines(vpds, rws["5%",  ], col = "#C75B39", lty = 2, lwd = 1.2)
lines(vpds, rws["50%", ], col = "#C75B39", lwd = 2.5)
