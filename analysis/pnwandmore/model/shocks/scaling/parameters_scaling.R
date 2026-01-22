rm(list = ls());gc()
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
library(terra)
library(rnaturalearth)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)


data_19001940 <-readRDS(file.path(wd, 'output/model', 'data_15jan2025_PIPO_standclimate_19001940.rds'))
fit_19001940 <- readRDS(file.path(wd, 'output/model/scaling', 'fit_15jan2025_PIPO_standclimate_19001940_customgrainsize.rds'))
samples_19001940 <- fit_19001940$draws(variables = c("beta_gdd", "beta_vpd", "beta_pre", "phi_sck", "rho_sp", "rho_sh"))

data_19401980 <- readRDS(file.path(wd, 'output/model', 'data_15jan2025_PIPO_standclimate_19401980.rds'))
fit_19401980 <- readRDS(file.path(wd, 'output/model/scaling', 'fit_15jan2025_PIPO_standclimate_19401980_customgrainsize.rds'))
samples_19401980 <- fit_19401980$draws(variables = c("beta_gdd", "beta_vpd", "beta_pre", "phi_sck", "rho_sp", "rho_sh"))

data_19802024 <- readRDS(file.path(wd, 'output/model', 'data_15jan2025_PIPO_standclimate_19802024.rds'))
fit_19802024 <- readRDS(file.path(wd, 'output/model/scaling', 'fit_15jan2025_PIPO_standclimate_19802024_customgrainsize.rds'))
samples_19802024 <- fit_19802024$draws(variables = c("beta_gdd", "beta_vpd", "beta_pre", "phi_sck", "rho_sp", "rho_sh"))

par(mfrow= c(4,1), mar = c(4,1,1,1), cex.lab = 1.4, cex.axis = 0.9)

util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,'beta_gdd[1]'], nrow = 1000, ncol = 4)), 100,
                                        flim = c(-0.15,0.15),
                                        display_name=expression(beta[gdd]),
                                        col = '#27278f')
util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,'beta_gdd[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[gdd]),
                                col = '#278f27', add = TRUE)
util$plot_expectand_pushforward(t(matrix(samples_19802024[1:1000,1:4,'beta_gdd[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[gdd]),
                                col = '#8f2727', add = TRUE)
legend(x = -0.15, y = 100,
       lwd = 2, 
       col = c( '#27278f', '#278f27', '#8f2727'),
       legend = c('1900-1940', '1940-1980', '1980-2024'),
       bty = 'n')


util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,'beta_vpd[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[vpd]),
                                col = '#27278f')
util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,'beta_vpd[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[vpd]),
                                col = '#278f27', add = TRUE)
util$plot_expectand_pushforward(t(matrix(samples_19802024[1:1000,1:4,'beta_vpd[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[vpd]),
                                col = '#8f2727', add = TRUE)

util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,'beta_pre[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[pre]),
                                col = '#27278f')
util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,'beta_pre[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[pre]),
                                col = '#278f27', add = TRUE)
util$plot_expectand_pushforward(t(matrix(samples_19802024[1:1000,1:4,'beta_pre[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(-0.15,0.15),
                                display_name=expression(beta[pre]),
                                col = '#8f2727', add = TRUE)


util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,'rho_sp[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(0,20),
                                display_name=expression(rho[species]),
                                col = '#27278f')
util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,'rho_sp[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(0,20),
                                display_name=expression(rho[species]),
                                col = '#278f27', add = TRUE)
util$plot_expectand_pushforward(t(matrix(samples_19802024[1:1000,1:4,'rho_sp[1]'], nrow = 1000, ncol = 4)), 100,
                                flim = c(0,20),
                                display_name=expression(rho[species]),
                                col = '#8f2727', add = TRUE)


util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,'rho_sh'], nrow = 1000, ncol = 4)), 100,
                                flim = c(0,20),
                                display_name=expression(rho[stand]),
                                col = '#27278f')
util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,'rho_sh'], nrow = 1000, ncol = 4)), 100,
                                flim = c(0,20),
                                display_name=expression(rho[stand]),
                                col = '#278f27', add = TRUE)




data_19001940$uniq_stand_ids %in% data_19401980$uniq_stand_ids

s <- 1
util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,paste0('phi_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,0.5),
                                display_name=expression(beta[pre]),
                                col = '#27278f')
util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,paste0('phi_sck[',which(data_19401980$uniq_stand_ids == data_19001940$uniq_stand_ids[s]),']')], nrow = 1000, ncol = 4)), 100,
                                flim = c(0,0.5),
                                display_name=expression(beta[pre]),
                                col = '#278f27', add = TRUE)


s <- 10
util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,paste0('phi_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,0.5),
                                display_name=expression(beta[pre]),
                                col = '#27278f')
util$plot_expectand_pushforward(t(matrix(samples_19401980[1:1000,1:4,paste0('phi_sck[',which(data_19401980$uniq_stand_ids == data_19001940$uniq_stand_ids[s]),']')], nrow = 1000, ncol = 4)), 100,
                                flim = c(0,0.5),
                                display_name=expression(beta[pre]),
                                col = '#278f27', add = TRUE)




util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,'phi_sck[1]'], nrow = 1000, ncol = 4)), 50,
                                flim = c(0,1),
                                display_name=expression(beta[gdd]),
                                col = '#27278f')
for(s in 2:data_19001940$N_stands){
  util$plot_expectand_pushforward(t(matrix(samples_19001940[1:1000,1:4,paste0('phi_sck[',s,']')], nrow = 1000, ncol = 4)), 50,
                                  flim = c(0,1),
                                  display_name=expression(beta[gdd]),
                                  col = '#27278f', add = TRUE)
}



stand_shocks <- data.frame(lat = data_19001940$uniq_stand_lat,
                           lon = data_19001940$uniq_stand_lon,
                           phi_sck = sapply(1:data_19001940$N_stands, function(s){
                             util$ensemble_mcmc_quantile_est(t(matrix(samples_19001940[1:1000,1:4,paste0('phi_sck[',s,']')])), c(0.5))}))

world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
westam <- vect(ne_states(country = c("Canada", "Mexico", "United States of America"), returnclass = "sf"))

plot(westam, xlim = c(-130,-100), ylim = c(30,50),
     col = 'grey90', border = 'white')
points(stand_shocks$lon, stand_shocks$lat, col = stand_shocks$phi_sck)


cols <- colorRampPalette(c("#5A9CB5", "#FAAC68", "#FA6868"))(length(seq(0,0.3,0.001)))

# map phi_sck to colors
phi_scaled <- cut(stand_shocks$phi_sck,
                  breaks = seq(0,0.3,0.001),
                  include.lowest = TRUE)

points(stand_shocks$lon,
       stand_shocks$lat,
       col = cols[phi_scaled],
       pch = 20)


## ---- color bar legend ----
z <- stand_shocks$phi_sck
par(xpd = TRUE)  # allow drawing outside plot

# position of color bar
x0 <- -129.5
x1 <- -128.5
y  <- seq(0,0.3,0.001)*20+30


rect(x0, y[-length(y)],
     x1, y[-1],
     col = cols,
     border = NA)

# axis + label
text(x = -129, y = 37, expression(phi[sck]))





# add numeric scale
z_ticks <- seq(0,0.3,0.1)
y_ticks <- z_ticks*20+30

axis(side = 4,
     at   = y_ticks,
     labels = z_ticks,
     las = 1,
     tck = -0.01,
     cex.axis = 0.8,
     pos = x1,  mgp = c(3, 0.5, 0))


par(xpd = FALSE)
