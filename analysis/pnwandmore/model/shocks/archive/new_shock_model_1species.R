rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)

data <- readRDS(file.path(wd, 'output/model', 'data_15nov2025_azca_pipo2018.rds'))
# fit <- readRDS(file.path(wd, "model/shocks/output", "fit_15nov2025_azca_pipo2018.rds"))
# samples <- util$extract_expectand_vals(fit)
samples <- readRDS(file.path(wd, "model/shocks/output", "samples_15nov2025_azca_pipo2018_gammash.rds"))

util$plot_expectand_pushforward(samples[['rho_sh']], 50,
                                flim = c(0,5), 
                                display_name=expression(rho[short]))

util$plot_expectand_pushforward(samples[['phi_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[shock]))

util$plot_expectand_pushforward(samples[['omega_tree_sck']], 50,
                                flim = c(0,1), 
                                display_name=expression(phi[shock]))


par(mfrow = c(3,2), mar = c(4,5,2,2), cex.main = 1)
random_trees <- sample(1:data$N_trees, 3)
for(t in random_trees){
  stand <- data$stand_idxs[t]
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  
  years_t <- data$all_years_idxs[idxs_t]
  names <- sapply(years_t,
                  function(x) paste0('f_sh[', stand,',',x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], 
    main = paste0('Stand ', stand), 
    ylab = expression(f[short]), display_ylim = c(-12,2), 
    display_xlim = data$all_years[c(1,length(data$all_years))])
  
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], 
    main = paste0('Stand ', stand, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), 
    ylab = expression(delta[shocks]), display_ylim = c(-12,2), 
    display_xlim = data$all_years[c(1,length(data$all_years))])
}  

nf <- layout(
  matrix(c(0,0,1,1,0,0,
           2,3,4,5,6,7,
           8,9,10,11,12,13), 
         ncol=6, byrow=TRUE)
)
stand <- 7
trees_stand <- which(data$stand_idxs == 7)
years <- 1:data$N_all_years
names <- sapply(years,
                function(x) paste0('f_sh[', stand,',',x, ']'))
util$plot_conn_pushforward_quantiles(
  samples, names, plot_xs = data$all_years, 
  main = paste0('Stand ', stand), 
  ylab = expression(f[short]), display_ylim = c(-6,2), 
  display_xlim = data$all_years[c(1,length(data$all_years))])
abline(v = 2002, lty = 2)
random_trees <- sample(trees_stand, 12)
for(t in random_trees){

  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  names <- sapply(idxs_t,
                  function(x) paste0('delta_sck[', x, ']'))
  util$plot_conn_pushforward_quantiles(
    samples, names, plot_xs = data$years[idxs_t], 
    main = paste0('Stand ', stand, ' - ', data$uniq_species_ids[data$species_idxs[t]], ' - tree ', t), 
    ylab = expression(delta[shocks]), display_ylim = c(-6,2), 
    display_xlim = data$all_years[c(1,length(data$all_years))])
  abline(v = 2002, lty = 2)
}  

# Check 2002
par(mfrow = c(3,3))
fsh <- paste0('f_sh[', stand,',',which(data$all_years == 2002), ']')
for(t in random_trees){
  idxs_t <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
  deltashock <- paste0('delta_sck[', idxs_t[data$years[idxs_t] == 2002], ']') 
  util$plot_pairs_by_chain(samples[[fsh]], expression(f[short]),
                           samples[[deltashock]], expression(delta[shock]))
}  

util$plot_pairs_by_chain(samples[['rho_sh']], expression(f[short]),
                         samples[['outer_tau_sck']], expression(delta[shock]))

