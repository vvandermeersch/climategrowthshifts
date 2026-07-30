par(mfrow = c(2,2))

# Northern Mountains
sites_here <- terra::intersect(sites, normnt)$id
all_years <- max(data$N_stand_years[sites_here])

n_sites_peryear <- rep(0, all_years)
for(i in 1:all_years){
  n_sites <- 0
  for(s in sites_here){
    shock <- gq_samples[[paste0('shock_growth_change[', s,',',i,']')]]
    if(all(shock == 'NaN')){
      next
    }else{
      n_sites <- n_sites + 1
    }
  }
  n_sites_peryear[i] <- n_sites
}

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(0,200), ylab= 'Number of sites',
     bty = 'n')
plot_xs <- data$all_years[1:all_years]
lines(plot_xs, n_sites_peryear, col="#7C873E", lwd=2)


# Sierra Nevada
sites_here <- terra::intersect(sites, sienev)$id
all_years <- max(data$N_stand_years[sites_here])

n_sites_peryear <- rep(0, all_years)
for(i in 1:all_years){
  n_sites <- 0
  for(s in sites_here){
    shock <- gq_samples[[paste0('shock_growth_change[', s,',',i,']')]]
    if(all(shock == 'NaN')){
      next
    }else{
      n_sites <- n_sites + 1
    }
  }
  n_sites_peryear[i] <- n_sites
}

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(0,200), ylab= 'Number of sites',
     bty = 'n')
plot_xs <- data$all_years[1:all_years]
lines(plot_xs, n_sites_peryear, col="#5495CF", lwd=2)

# Colorado plateau and more
sites_here <- terra::intersect(sites, colplat)$id
all_years <- max(data$N_stand_years[sites_here])

n_sites_peryear <- rep(0, all_years)
for(i in 1:all_years){
  n_sites <- 0
  for(s in sites_here){
    shock <- gq_samples[[paste0('shock_growth_change[', s,',',i,']')]]
    if(all(shock == 'NaN')){
      next
    }else{
      n_sites <- n_sites + 1
    }
  }
  n_sites_peryear[i] <- n_sites
}

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(0,200), ylab= 'Number of sites',
     bty = 'n')
plot_xs <- data$all_years[1:all_years]
lines(plot_xs, n_sites_peryear, col="#DB4743", lwd=2)

# Deserts
sites_here <- terra::intersect(sites, deserts)$id
all_years <- max(data$N_stand_years[sites_here])

n_sites_peryear <- rep(0, all_years)
for(i in 1:all_years){
  n_sites <- 0
  for(s in sites_here){
    shock <- gq_samples[[paste0('shock_growth_change[', s,',',i,']')]]
    if(all(shock == 'NaN')){
      next
    }else{
      n_sites <- n_sites + 1
    }
  }
  n_sites_peryear[i] <- n_sites
}

plot(1, type="n", main='',
     xlim=range(data$all_years), xlab='',
     ylim=c(0,200), ylab= 'Number of sites',
     bty = 'n')
plot_xs <- data$all_years[1:all_years]
lines(plot_xs, n_sites_peryear, col="#F5AF4D", lwd=2)

