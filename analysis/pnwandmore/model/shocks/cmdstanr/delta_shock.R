library(cmdstanr)
model <- cmdstan_model(file.path(wd, 'model/stan', 'model10_speciesshocks_improved_nopooling_standclimate_fshinthreads_grainsize_wgenquant.stan'),
                       cpp_options = list(stan_threads = TRUE))

data <- data_19001940


data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
data$uniq_stand_lon <- NULL
data$uniq_stand_lat <- NULL
data$N_clades <- 1

data$clade_idxs <- array(data$clade_idxs, dim = 1)

data$grainsize <- ceiling(data$N_stands/16)



fit <- model$sample(data = data, threads_per_chain = 16,
                    chains = 1, parallel_chains = 1, seed = 23079696,
                    refresh = 5, iter_warmup = 10, iter_sampling = 100)

delta_sck <- fit$draws(variables = c("delta_sck"))

t <- 3
tidxs <- data$tree_start_idxs[t]:data$tree_end_idxs[t]
delta_sck_t <- delta_sck[,,paste0('delta_sck[',tidxs,']')]
delta_sck_t <- lapply(1:dim(delta_sck_t)[3], function(i){t(matrix(delta_sck_t[,,i], nrow = dim(delta_sck_t)[1], ncol = dim(delta_sck_t)[2]))})
names(delta_sck_t) <- paste0('delta_sck[',tidxs,']')


util$plot_conn_pushforward_quantiles(
  delta_sck_t, names(delta_sck_t), plot_xs = data$years[tidxs], 
  main = '', 
  ylab = expression(delta[shocks]), display_ylim = c(-5,5), 
  display_xlim = data$all_years[c(1,length(data$all_years))])
