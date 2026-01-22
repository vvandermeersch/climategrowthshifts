


wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
# setwd(wd)
# util <- new.env()
# source('mcmc_analysis_tools_rstan.R', local=util)
# source('mcmc_visualization_tools.R', local=util)


library(cmdstanr)

data <- readRDS(file.path(wd, 'output/model', 'data_28dec2025_8standsPIPO_standclimate.rds'))
data$uniq_tree_ids <- NULL
data$uniq_species_ids <- NULL
data$uniq_stand_ids <- NULL
data$N_clades <- 1

data$clade_idxs <- array(data$clade_idxs, dim = 1)

model <- cmdstan_model(file.path(wd, 'model/stan', 'model10_speciesshocks_improved_nopooling_standclimate_fshinthreads.stan'),
                       cpp_options = list(stan_threads = TRUE))

fit <- model$sample(data = data, threads_per_chain = 2,
                    chains = 4, parallel_chains = 4, seed = 23079696,
                    refresh = 10, iter_warmup = 50, iter_sampling = 50)

profiling <- fit$profiles()
runtime <- fit$time()

fit2 <- model$sample(data = data, threads_per_chain = 4,
                    chains = 4, parallel_chains = 4, seed = 23079696,
                    refresh = 10, iter_warmup = 50, iter_sampling = 50)

profiling2 <- fit2$profiles()
runtime2 <- fit2$time()


plot(x = NULL,
     y = NULL,
     xlim = c(0, 54),
     ylim = c(0, 30*60),
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = 'Runtime (minutes)',
     bty = "n")


axis(1, at = seq(2.5,20,5), labels = paste0('Chain ', 1:4))
axis(2, at = seq(0, 30*60, 5*60), labels = seq(0, 30*60, 5*60)/60)

for(c in 1:length(profiling)){
  
  profile_chain <- data.frame(profiling[[c]])
  threads <- unique(profile_chain$thread_id)
  
  for(t in 1:length(threads)){
    
    profile_thread <- profile_chain[profile_chain$thread_id == threads[t],]
    
    rect(xleft = 0.5+(t-1)+(c-1)*5,
         ybottom = 0,
         xright = 0.5+t+(c-1)*5,
         ytop =  profile_thread[profile_thread$name == 'compute_f', 'total_time'],
         col =  '#F6CF71',
         border = 'white') 
    sumtime <- profile_thread[profile_thread$name == 'compute_f', 'total_time']
    
    rect(xleft = 0.5+(t-1)+(c-1)*5,
         ybottom = sumtime,
         xright = 0.5+t+(c-1)*5,
         ytop =  sumtime + profile_thread[profile_thread$name == 'compute_mu', 'total_time'],
         col =  '#F89C74',
         border = 'white') 
    sumtime <- profile_thread[profile_thread$name == 'compute_f', 'total_time'] + 
      profile_thread[profile_thread$name == 'compute_mu', 'total_time']
    
    rect(xleft = 0.5+(t-1)+(c-1)*5,
         ybottom = sumtime,
         xright = 0.5+t+(c-1)*5,
         ytop =  sumtime + profile_thread[profile_thread$name == 'store_trees', 'total_time'],
         col =  '#87C55F',
         border = 'white') 
    sumtime <- profile_thread[profile_thread$name == 'compute_f', 'total_time'] + 
      profile_thread[profile_thread$name == 'compute_mu', 'total_time'] +
      profile_thread[profile_thread$name == 'store_trees', 'total_time']
    
    rect(xleft = 0.5+(t-1)+(c-1)*5,
         ybottom = sumtime,
         xright = 0.5+t+(c-1)*5,
         ytop =  sumtime + profile_thread[profile_thread$name == 'compute_fsh', 'total_time'],
         col =  '#11A579',
         border = 'white') 

    rect(xleft = 0.5+(t-1)+(c-1)*5,
         ybottom = profile_thread[profile_thread$name == 'slice_loop', 'total_time'],
         xright = 0.5+t+(c-1)*5,,
         ytop =  profile_chain[profile_chain$name == 'likelihood', 'total_time'],
         col =  NA,
         density = 20,
         angle = 45,
         border = NA) 
    
    if(c == 1){
      text(x =  (t-1)+(c-1)*5+1, y = 30, labels = paste0('Thread ', t), srt = 90, adj = 0, cex = 0.8)
    }
    
  }
  
  
  rect(xleft = 0.5+(c-1)*5,
       ybottom = profile_chain[profile_chain$name == 'likelihood', 'total_time'],
       xright = 0.5+length(threads)+(c-1)*5,
       ytop =  profile_chain[profile_chain$name == 'likelihood', 'total_time'] +
         profile_chain[profile_chain$name == 'L_cov_fsh', 'total_time'],
       col =  '#DCB0F2',
       border = NA)
  
  rect(xleft = 0.5+(c-1)*5,
       ybottom = profile_chain[profile_chain$name == 'likelihood', 'total_time'] +
         profile_chain[profile_chain$name == 'L_cov_fsh', 'total_time'],
       xright = 0.5+length(threads)+(c-1)*5,
       ytop =  profile_chain[profile_chain$name == 'likelihood', 'total_time'] +
         profile_chain[profile_chain$name == 'L_cov_fsh', 'total_time'] +
         profile_chain[profile_chain$name == 'f_lg', 'total_time'],
       col =  '#9EB9F3',
       border = NA)
  
  rect(xleft = 0.5+(c-1)*5,
       ybottom = 0,
       xright = 0.5+length(threads)+(c-1)*5,
       ytop =  runtime$chains[c, 'total'],
       col =  NA,
       border = 'black') 
  
}

# segments(x0 = 1+(4-1)+(4-1)*5, x1 = 23, y0 = 17*60, y1 = 17*60)
# text(x = 23.5, y = 17*60, labels = 'Overhead?', srt = 0, adj = 0, cex = 0.8)
# 
# segments(x0 = 1+(4-1)+(4-1)*5, x1 = 23, y0 = 12*60, y1 = 12*60)
# text(x = 23.5, y = 12*60, labels = 'Loop over years and\nstore tree-level shocks', srt = 0, adj = 0, cex = 0.8)
# 
# segments(x0 = 1+(4-1)+(4-1)*5, x1 = 23, y0 = 8*60, y1 = 8*60)
# text(x = 23.5, y = 8*60, labels = 'Sum of climate predictors', srt = 0, adj = 0, cex = 0.8)
# 
# segments(x0 = 1+(4-1)+(4-1)*5, x1 = 23, y0 = 4*60, y1 = 4*60)
# text(x = 23.5, y = 4*60, labels = 'Compute tree-level GP', srt = 0, adj = 0, cex = 0.8)
# 
# segments(x0 = 1+(4-1)+(4-1)*5, x1 = 23, y0 = 23*60, y1 = 23*60)
# text(x = 23.5, y = 23*60, labels = 'Compute stand-level GP', srt = 0, adj = 0, cex = 0.8)
# 
# segments(x0 = 1+(4-1)+(4-1)*5, x1 = 23, y0 = 25.5*60, y1 = 25.5*60)
# text(x = 23.5, y = 25.5*60, labels = 'Cholesky decomposition\nfor tree-level GP', srt = 0, adj = 0, cex = 0.8)


# plot(x = NULL,
#      y = NULL,
#      xlim = c(0, 50),
#      ylim = c(0, 35*60),
#      xaxt = 'n',
#      yaxt = 'n',
#      xlab = '',
#      ylab = 'Runtime (minutes)',
#      bty = "n")


axis(1, at = 20+seq(2.5,20,5), labels = paste0('Chain ', 1:4))
#axis(2, at = seq(0, 30*60, 5*60), labels = seq(0, 30*60, 5*60)/60)

for(c in 1:length(profiling)){
  
  profile_chain <- data.frame(profiling2[[c]])
  threads <- unique(profile_chain$thread_id)
  
  for(t in 1:length(threads)){
    
    profile_thread <- profile_chain[profile_chain$thread_id == threads[t],]
    
    rect(xleft = 20+0.5+(t-1)+(c-1)*5,
         ybottom = 0,
         xright =  20+0.5+t+(c-1)*5,
         ytop =  profile_thread[profile_thread$name == 'compute_f', 'total_time'],
         col =  '#F6CF71',
         border = 'white') 
    sumtime <- profile_thread[profile_thread$name == 'compute_f', 'total_time']
    
    rect(xleft =  20+0.5+(t-1)+(c-1)*5,
         ybottom = sumtime,
         xright =  20+0.5+t+(c-1)*5,
         ytop =  sumtime + profile_thread[profile_thread$name == 'compute_mu', 'total_time'],
         col =  '#F89C74',
         border = 'white') 
    sumtime <- profile_thread[profile_thread$name == 'compute_f', 'total_time'] + 
      profile_thread[profile_thread$name == 'compute_mu', 'total_time']
    
    rect(xleft =  20+0.5+(t-1)+(c-1)*5,
         ybottom = sumtime,
         xright =  20+0.5+t+(c-1)*5,
         ytop =  sumtime + profile_thread[profile_thread$name == 'store_trees', 'total_time'],
         col =  '#87C55F',
         border = 'white') 
    sumtime <- profile_thread[profile_thread$name == 'compute_f', 'total_time'] + 
      profile_thread[profile_thread$name == 'compute_mu', 'total_time'] +
      profile_thread[profile_thread$name == 'store_trees', 'total_time']
    
    rect(xleft = 20+0.5+(t-1)+(c-1)*5,
         ybottom = sumtime,
         xright = 20+0.5+t+(c-1)*5,
         ytop =  sumtime + profile_thread[profile_thread$name == 'compute_fsh', 'total_time'],
         col =  '#11A579',
         border = 'white') 
    
    rect(xleft =  20+0.5+(t-1)+(c-1)*5,
         ybottom = profile_thread[profile_thread$name == 'slice_loop', 'total_time'],
         xright =  20+0.5+t+(c-1)*5,,
         ytop =  profile_chain[profile_chain$name == 'likelihood', 'total_time'],
         col =  NA,
         density = 20,
         angle = 45,
         border = NA) 
    
    if(c == 1){
      text(x =  20+(t-1)+(c-1)*5+1, y = 30, labels = paste0('Thread ', t), srt = 90, adj = 0, cex = 0.8)
    }
    
  }
  
  
  rect(xleft =  20+0.5+(c-1)*5,
       ybottom = profile_chain[profile_chain$name == 'likelihood', 'total_time'],
       xright =  20+0.5+length(threads)+(c-1)*5,
       ytop =  profile_chain[profile_chain$name == 'likelihood', 'total_time'] +
         profile_chain[profile_chain$name == 'L_cov_fsh', 'total_time'],
       col =  '#DCB0F2',
       border = NA)
  
  rect(xleft =  20+0.5+(c-1)*5,
       ybottom = profile_chain[profile_chain$name == 'likelihood', 'total_time'] +
         profile_chain[profile_chain$name == 'L_cov_fsh', 'total_time'],
       xright =  20+0.5+length(threads)+(c-1)*5,
       ytop =  profile_chain[profile_chain$name == 'likelihood', 'total_time'] +
         profile_chain[profile_chain$name == 'L_cov_fsh', 'total_time'] +
         profile_chain[profile_chain$name == 'f_lg', 'total_time'],
       col =  '#9EB9F3',
       border = NA)
  
  rect(xleft =  20+0.5+(c-1)*5,
       ybottom = 0,
       xright =  20+0.5+length(threads)+(c-1)*5,
       ytop =  runtime2$chains[c, 'total'],
       col =  NA,
       border = 'black') 
  
}


segments(x0 = 20+1+(4-1)+(4-1)*5, x1 =  20+23, y0 = 19*60, y1 = 19*60)
text(x =  20+23.5, y = 19*60, labels = 'Overhead?', srt = 0, adj = 0, cex = 0.8)

segments(x0 =  20+1+(4-1)+(4-1)*5, x1 =  20+23, y0 = 15.1*60, y1 = 15.1*60)
text(x =  20+23.5, y = 15.1*60, labels = 'Compute stand-level GP', srt = 0, adj = 0, cex = 0.8)

segments(x0 =  20+1+(4-1)+(4-1)*5, x1 =  20+23, y0 = 11*60, y1 = 11*60)
text(x =  20+23.5, y = 11*60, labels = 'Loop over years and\nstore tree-level shocks', srt = 0, adj = 0, cex = 0.8)

segments(x0 =  20+1+(4-1)+(4-1)*5, x1 =  20+23, y0 = 8*60, y1 = 8*60)
text(x =  20+23.5, y = 8*60, labels = 'Sum of climate predictors', srt = 0, adj = 0, cex = 0.8)

segments(x0 =  20+1+(4-1)+(4-1)*5, x1 =  20+23, y0 = 4*60, y1 = 4*60)
text(x =  20+23.5, y = 4*60, labels = 'Compute tree-level GP', srt = 0, adj = 0, cex = 0.8)


segments(x0 =  20+1+(4-1)+(4-1)*5, x1 =  20+23, y0 = 25.5*60, y1 = 25.5*60)
text(x =  20+23.5, y = 26.5*60, labels = 'Cholesky decomposition\nfor tree-level GP', srt = 0, adj = 0, cex = 0.8)

segments(x0 =  20+1+(4-1)+(4-1)*5, x1 =  20+23, y0 = 23.7*60, y1 = 23.7*60)
text(x =  20+23.5, y = 23.7*60, labels = 'for stand-level GP', srt = 0, adj = 0, cex = 0.8)

text(x =5, y = 30*60, labels = '2 threads per chain', srt = 0, adj = 0, cex = 0.8)
text(x =25, y = 30*60, labels = '4 threads per chain', srt = 0, adj = 0, cex = 0.8)
