rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_11july2025_soilmoisture2m.rds'))

# Extract predictions for one species
speciesid <- 'PSME'
focustrees <- which(data$species_idxs == which(data$uniq_species_ids == speciesid))
focuspred <- unlist(sapply(focustrees, function(t)
  return(data$tree_start_idxs[t]:data$tree_end_idxs[t])))
prednames <- sapply(focuspred, function(n) paste0('log_rw_pred[', n, ']'))
full <- FALSE
if(full){
  fit <- readRDS(file.path(wd, "output/model/fit_11july2025_partialpooling_2clades_centered_extended_soilmoisture2m.rds"))
  
  # samples <- util$extract_expectand_vals(fit)
  predictions <- rstan::extract(fit, pars = prednames, permuted=FALSE)
  # Below, from Mike's function
  N <- length(prednames)
  predictions <- lapply(1:N, function(n) t(predictions[,,n]))
  names(predictions) <- prednames
  
  saveRDS(predictions, file.path(wd, "output/model", paste0(speciesid, "_predictions_11july2025_partialpooling_2clades_centered_extended_soilmoisture2m.rds")))
  
}else{
  predictions <- readRDS(file.path(wd, "output/model", paste0(speciesid, "_predictions_11july2025_partialpooling_2clades_centered_extended_soilmoisture2m.rds")))
}

# Let's take 4 first PIPO trees (all one the same stand)
data$uniq_stand_ids[data$stand_idxs[focustrees[2601:2604]]]
for(t in focustrees[2701:2704]){
  print(t)
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  prednames <- sapply(idxs, function(n) paste0('log_rw_pred[', n, ']'))
  util$plot_conn_pushforward_quantiles(predictions, prednames, data$years[idxs],
                                       xlab="Year", ylab="Log Ring Width Per mm", 
                                       display_ylim=c(-3, 3), display_xlim = c(1980, 2024))
  # points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.0, col="white")
  # points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col="black")
}

# The predictions are very similar
par(mfrow = c(1,1))
list_colors <- list(
  green_colors <- c("#A0CFA0", "#76B876", "#66A766", 
                    "#4F9F4F", "#388E38", "#2B7A2B"),
  teal_colors <- c("#A3D9D0", "#80B7B0", "#6E9C9C", 
                   "#5A8080", "#4A6D6D", "#395757"),
  blue_colors <- c("#A0C0D9", "#7A98B8", "#6886A6", 
                   "#4F6D94", "#38608A", "#2A4B73"),
  purple_colors <- c("#D4A0D9", "#B97AB0", "#A66E9F", 
                     "#8F528E", "#7A2D7A", "#6A0066")
)
col <- 1
N_plots = 10;plot_xs =data$years[idxs]
display_xlim = c(1980, 2010)
display_ylim=c(-3, 2)
xlab="Year";ylab="Log Ring Width Per mm"
par(mar = c(4,4,1,1))
plot(1, type="n", main='',
     xlim=display_xlim, xlab=xlab,
     ylim=display_ylim, ylab=ylab)
for(t in focustrees[2701:2704]){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  prednames <- sapply(idxs, function(n) paste0('log_rw_pred[', n, ']'))
  
  fs <- t(sapply(prednames, function(name)
    c(predictions[[name]], recursive=TRUE)))
  N <- dim(fs)[1]
  I <- dim(fs)[2]
  
  # Configure ensemble of function realizations
  J <- min(N_plots, I)
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  line_colors <- colormap(colormap=list_colors[[col]], nshades=J)
  
  # Plot
  delta <- 0.05 * (display_ylim[2] - display_ylim[1])
  display_ylim[1] <- display_ylim[1] - delta
  display_ylim[2] <- display_ylim[2] + delta
  
  
  
  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    lines(plot_xs, f, col=line_colors[j], lwd=1)
  }

  col <- col + 1
}
col <- 1
for(t in focustrees[2701:2704]){
  idxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=1.2, col="white")
  points(data$years[idxs], data$log_rw_obs[idxs], pch=16, cex=0.8, col= list_colors[[col]][6])
  
  col <- col + 1
}
abline(v=2005, lwd=1, lty=2, col="grey50")

# Look at climate on this stand
par(mfrow = c(3,1), mar = c(4,4,1,1))
display_ylim=c(5, 17)
plot(1, type="n", main='',
     xlim=display_xlim, xlab=xlab,
     ylim=display_ylim, ylab='VPD (hPa)')
lines(data$years[idxs], data$vpd_obs[idxs], pch=16, cex=0.8, col= "black")
points(data$years[idxs], data$vpd_obs[idxs], pch=16, cex=1.2, col="white")
points(data$years[idxs], data$vpd_obs[idxs], pch=16, cex=0.8, col= "black")
abline(v=2005, lwd=1, lty=2, col="grey50")


display_ylim=c(10, 22)
plot(1, type="n", main='',
     xlim=display_xlim, xlab=xlab,
     ylim=display_ylim, ylab='GDD (x100 degC)')
lines(data$years[idxs], data$gdd_obs[idxs], pch=16, cex=0.8, col= "black")
points(data$years[idxs], data$gdd_obs[idxs], pch=16, cex=1.2, col="white")
points(data$years[idxs], data$gdd_obs[idxs], pch=16, cex=0.8, col= "black")
abline(v=2005, lwd=1, lty=2, col="grey50")

display_ylim=c(18, 28)
plot(1, type="n", main='',
     xlim=display_xlim, xlab=xlab,
     ylim=display_ylim, ylab='SM (%)')
lines(data$years[idxs], data$sm_obs[idxs], pch=16, cex=0.8, col= "black")
points(data$years[idxs], data$sm_obs[idxs], pch=16, cex=1.2, col="white")
points(data$years[idxs], data$sm_obs[idxs], pch=16, cex=0.8, col= "black")
abline(v=2005, lwd=1, lty=2, col="grey50")

par(mfrow = c(1,1), mar = c(4,4,1,1))
stand <- unique(data$stand_idxs[focustrees[20:23]])
prednames <- sapply(1:data$N_all_years, function(n) paste0('f_sh[', stand,',', n, ']'))
samples_fsh <- rstan::extract(fit, pars = prednames, permuted=FALSE)
N <- length(prednames)
samples_fsh <- lapply(1:N, function(n) t(samples_fsh[,,n]))
names(samples_fsh) <- prednames
util$plot_realizations(samples_fsh, prednames, data$all_years, N_plots=50,
                       xlab="Year", display_ylim=c(-3, 3), display_xlim = c(1980, 2001))
abline(v=1996, lwd=1, lty=2, col="grey50")
abline(v=1992, lwd=1, lty=2, col="grey50")
  
