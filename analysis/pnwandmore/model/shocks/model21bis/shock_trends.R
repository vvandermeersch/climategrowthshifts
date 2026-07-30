par(mfrow = c(2,2))

# Northern Mountains
sites_here <- terra::intersect(sites, normnt)$id
for(i in 1:data$N_all_years){
  base_samples[[paste0('mean_phi[',i,']')]] <- 0
}
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[[paste0('mean_phi[',i,']')]] <- base_samples[[paste0('mean_phi[',i,']')]] + phi/length(sites_here)
  }
}

names <- paste0('mean_phi[',1:data$N_all_years,']')
plot_xs <- data$all_years
# Construct quantiles for bin contents
N <- length(names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(base_samples[[names[n]]], probs)
}
plot_quantiles <- sapply(1:N, calc)

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(0, 0.5), ylab= 'Avg. p(concordant event)',
     btry = 'n')

polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[1,], rev(plot_quantiles[9,])),
        col = "#7C873E40", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[2,], rev(plot_quantiles[8,])),
        col = "#7C873E60", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[3,], rev(plot_quantiles[7,])),
        col = "#7C873E80", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[4,], rev(plot_quantiles[6,])),
        col = "#7C873E80", border = NA)
lines(plot_xs, plot_quantiles[5,], col="#7C873E", lwd=2)
before2000 <- mean(sapply(which(data$all_years < 2000), 
                          function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
after2000 <- mean(sapply(which(data$all_years >= 2000), 
                         function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
segments(x0 = 1950, x1 = 1999, y0 = before2000, lty = 2, col = "#7C873E")
segments(x0 = 2000, x1 = 2024, y0 = after2000, lty = 2, col = "#7C873E")

# Sierra Nevada
sites_here <- terra::intersect(sites, sienev)$id
for(i in 1:data$N_all_years){
  base_samples[[paste0('mean_phi[',i,']')]] <- 0
}
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[[paste0('mean_phi[',i,']')]] <- base_samples[[paste0('mean_phi[',i,']')]] + phi/length(sites_here)
  }
}

names <- paste0('mean_phi[',1:data$N_all_years,']')
plot_xs <- data$all_years
# Construct quantiles for bin contents
N <- length(names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(base_samples[[names[n]]], probs)
}
plot_quantiles <- sapply(1:N, calc)

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(0, 0.5), ylab= 'Avg. p(concordant event)',
     btry = 'n')

polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[1,], rev(plot_quantiles[9,])),
        col = "#5495CF40", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[2,], rev(plot_quantiles[8,])),
        col = "#5495CF60", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[3,], rev(plot_quantiles[7,])),
        col = "#5495CF80", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[4,], rev(plot_quantiles[6,])),
        col = "#5495CF80", border = NA)
lines(plot_xs, plot_quantiles[5,], col="#5495CF", lwd=2)
before2000 <- mean(sapply(which(data$all_years < 2000), 
                          function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
after2000 <- mean(sapply(which(data$all_years >= 2000), 
                         function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
segments(x0 = 1950, x1 = 1999, y0 = before2000, lty = 2, col = "#5495CF")
segments(x0 = 2000, x1 = 2024, y0 = after2000, lty = 2, col = "#5495CF")

# Colorado plateau and more
sites_here <- terra::intersect(sites, colplat)$id
for(i in 1:data$N_all_years){
  base_samples[[paste0('mean_phi[',i,']')]] <- 0
}
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[[paste0('mean_phi[',i,']')]] <- base_samples[[paste0('mean_phi[',i,']')]] + phi/length(sites_here)
  }
}

names <- paste0('mean_phi[',1:data$N_all_years,']')
plot_xs <- data$all_years
# Construct quantiles for bin contents
N <- length(names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(base_samples[[names[n]]], probs)
}
plot_quantiles <- sapply(1:N, calc)

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(0, 0.5), ylab= 'Avg. p(concordant event)',
     btry = 'n')

polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[1,], rev(plot_quantiles[9,])),
        col = "#DB474340", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[2,], rev(plot_quantiles[8,])),
        col = "#DB474360", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[3,], rev(plot_quantiles[7,])),
        col = "#DB474380", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[4,], rev(plot_quantiles[6,])),
        col = "#DB474380", border = NA)
lines(plot_xs, plot_quantiles[5,], col="#DB4743", lwd=2)
before2000 <- mean(sapply(which(data$all_years < 2000), 
                          function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
after2000 <- mean(sapply(which(data$all_years >= 2000), 
                         function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
segments(x0 = 1950, x1 = 1999, y0 = before2000, lty = 2, col = "#DB4743")
segments(x0 = 2000, x1 = 2024, y0 = after2000, lty = 2, col = "#DB4743")

# Deserts
sites_here <- terra::intersect(sites, deserts)$id
for(i in 1:data$N_all_years){
  base_samples[[paste0('mean_phi[',i,']')]] <- 0
}
for(s in sites_here){
  years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)
  for(i in 1:length(years)){
    
    logit_phi <- base_samples[[paste0('logit_phi_sck[', s, ']')]]+
      (data$vpd_obs[years[i]]-vpd0)*base_samples[['beta_phi_vpd']]+
      (data$pre_obs[years[i]]-pre0)*base_samples[['beta_phi_pre']]
    phi <- boot::inv.logit(logit_phi)
    
    base_samples[[paste0('mean_phi[',i,']')]] <- base_samples[[paste0('mean_phi[',i,']')]] + phi/length(sites_here)
  }
}

names <- paste0('mean_phi[',1:data$N_all_years,']')
plot_xs <- data$all_years
# Construct quantiles for bin contents
N <- length(names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(base_samples[[names[n]]], probs)
}
plot_quantiles <- sapply(1:N, calc)

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(0, 0.5), ylab= 'Avg. p(concordant event)',
     btry = 'n')

polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[1,], rev(plot_quantiles[9,])),
        col = "#F5AF4D40", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[2,], rev(plot_quantiles[8,])),
        col = "#F5AF4D60", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[3,], rev(plot_quantiles[7,])),
        col = "#F5AF4D80", border = NA)
polygon(c(plot_xs, rev(plot_xs)),
        c(plot_quantiles[4,], rev(plot_quantiles[6,])),
        col = "#F5AF4D80", border = NA)
lines(plot_xs, plot_quantiles[5,], col="#F5AF4D", lwd=2)
before2000 <- mean(sapply(which(data$all_years < 2000), 
                          function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
after2000 <- mean(sapply(which(data$all_years >= 2000), 
                         function(y) util$ensemble_mcmc_quantile_est(base_samples[[paste0('mean_phi[',y,']')]], c(0.5))))
segments(x0 = 1950, x1 = 1999, y0 = before2000, lty = 2, col = "#F5AF4D")
segments(x0 = 2000, x1 = 2024, y0 = after2000, lty = 2, col = "#F5AF4D")
