par(mfrow = c(2,2))

# Northern Mountains
sites_here <- terra::intersect(sites, normnt)$id
all_years <- max(data$N_stand_years[sites_here])
for(i in 1:all_years){
  gq_samples[[paste0('mean_shock[',i,']')]] <- 0
}

for(i in 1:all_years){
  n_sites <- 0
  for(s in sites_here){
    shock <- gq_samples[[paste0('shock_growth_change[', s,',',i,']')]]
    if(all(shock == 'NaN')){
      next
    }else{
      gq_samples[[paste0('mean_shock[',i,']')]] <- gq_samples[[paste0('mean_shock[',i,']')]] + shock
      n_sites <- n_sites + 1
    }
  }
  gq_samples[[paste0('mean_shock[',i,']')]] <- gq_samples[[paste0('mean_shock[',i,']')]]/n_sites
}

names <- paste0('mean_shock[',1:all_years,']')
plot_xs <- data$all_years[1:all_years]
# Construct quantiles for bin contents
N <- length(names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(gq_samples[[names[n]]], probs)
}
plot_quantiles <- sapply(1:N, calc)

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(-45, 5), ylab= 'Average growth variation\ndue to concordant shocks (%)',
     bty = 'n')

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


# Sierra Nevada
sites_here <- terra::intersect(sites, sienev)$id
all_years <- max(data$N_stand_years[sites_here])
for(i in 1:all_years){
  gq_samples[[paste0('mean_shock[',i,']')]] <- 0
}

for(i in 1:all_years){
  n_sites <- 0
  for(s in sites_here){
    shock <- gq_samples[[paste0('shock_growth_change[', s,',',i,']')]]
    if(all(shock == 'NaN')){
      next
    }else{
      gq_samples[[paste0('mean_shock[',i,']')]] <- gq_samples[[paste0('mean_shock[',i,']')]] + shock
      n_sites <- n_sites + 1
    }
  }
  gq_samples[[paste0('mean_shock[',i,']')]] <- gq_samples[[paste0('mean_shock[',i,']')]]/n_sites
}

names <- paste0('mean_shock[',1:all_years,']')
plot_xs <- data$all_years[1:all_years]
# Construct quantiles for bin contents
N <- length(names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(gq_samples[[names[n]]], probs)
}
plot_quantiles <- sapply(1:N, calc)

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(-45, 5), ylab= 'Average growth variation\ndue to concordant shocks (%)',
     bty = 'n')

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

# Colorado plateau and more
sites_here <- terra::intersect(sites, colplat)$id
all_years <- max(data$N_stand_years[sites_here])
for(i in 1:all_years){
  gq_samples[[paste0('mean_shock[',i,']')]] <- 0
}

for(i in 1:all_years){
  n_sites <- 0
  for(s in sites_here){
    shock <- gq_samples[[paste0('shock_growth_change[', s,',',i,']')]]
    if(all(shock == 'NaN')){
      next
    }else{
      gq_samples[[paste0('mean_shock[',i,']')]] <- gq_samples[[paste0('mean_shock[',i,']')]] + shock
      n_sites <- n_sites + 1
    }
  }
  gq_samples[[paste0('mean_shock[',i,']')]] <- gq_samples[[paste0('mean_shock[',i,']')]]/n_sites
}

names <- paste0('mean_shock[',1:all_years,']')
plot_xs <- data$all_years[1:all_years]
# Construct quantiles for bin contents
N <- length(names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(gq_samples[[names[n]]], probs)
}
plot_quantiles <- sapply(1:N, calc)

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(-45, 5), ylab= 'Average growth variation\ndue to concordant shocks (%)',
     bty = 'n')

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

# Deserts
sites_here <- terra::intersect(sites, deserts)$id
all_years <- max(data$N_stand_years[sites_here])
for(i in 1:all_years){
  gq_samples[[paste0('mean_shock[',i,']')]] <- 0
}

for(i in 1:all_years){
  n_sites <- 0
  for(s in sites_here){
    shock <- gq_samples[[paste0('shock_growth_change[', s,',',i,']')]]
    if(all(shock == 'NaN')){
      next
    }else{
      gq_samples[[paste0('mean_shock[',i,']')]] <- gq_samples[[paste0('mean_shock[',i,']')]] + shock
      n_sites <- n_sites + 1
    }
  }
  gq_samples[[paste0('mean_shock[',i,']')]] <- gq_samples[[paste0('mean_shock[',i,']')]]/n_sites
}

names <- paste0('mean_shock[',1:all_years,']')
plot_xs <- data$all_years[1:all_years]
# Construct quantiles for bin contents
N <- length(names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(gq_samples[[names[n]]], probs)
}
plot_quantiles <- sapply(1:N, calc)

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(-45, 5), ylab= 'Average growth variation\ndue to concordant shocks (%)',
     bty = 'n')

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

