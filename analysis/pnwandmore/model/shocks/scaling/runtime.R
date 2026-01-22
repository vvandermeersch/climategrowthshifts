wd <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore'
library(rstan)
library(cmdstanr)

fit8stands <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO.rds'))
fit8stands_cmdstand <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO_cmdstan.rds'))
# fit16stands <- readRDS(file.path(wd, 'output/model/scaling', 'fit_16standsPIPO.rds'))
fit8stands_ov <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO_oldversion.rds'))
fit8stands_standpred <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO_standpred.rds'))
fit8stands_standpred_improved <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO_standpred_improved.rds'))
fit8stands_standpred_improveds_cmdstand <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO_standpred_improved_cmdstan.rds'))

plot(x = NULL,
     y = NULL,
     xlim = c(0, 30),
     ylim = c(0, 3600*4),
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = 'Runtime (hours)',
     bty = "n")
par(mgp = c(3, 1.5, 0), mfrow = c(1,1))

axis(1, at = seq(2.5,27.5,5), labels = c('Initial\nversion', '4 cores\n(Rstan)', '4 cores\n(cmdstan)',
                                          'Stand pred.\n(Rstan)', 'Stand pred.\nv2 (Rstan)', 'Stand pred.\nv2 (cmdstan)'))
axis(2, at = seq(0, 3600*4, 3600), labels = sapply(seq(0, 4, 1), function(i) ifelse(i%%2==0,i,'')))

k <- 0
time <- get_elapsed_time(fit8stands_ov)
for(i in 1:4){
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = 0,
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white') 
  
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = time[i,'warmup'],
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup']+time[i,'sample'],
       col =  '#87C55F',
       border = 'white') 
}

k <- 1
time <- get_elapsed_time(fit8stands)
for(i in 1:4){
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = 0,
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white') 
  
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = time[i,'warmup'],
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup']+time[i,'sample'],
       col =  '#87C55F',
       border = 'white') 
}

k <- 2
time <- fit8stands_cmdstand$time()$chains
for(i in 1:4){
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = 0,
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white') 
  
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = time[i,'warmup'],
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup']+time[i,'sampling'],
       col =  '#87C55F',
       border = 'white') 
}

k <- 3
time <- get_elapsed_time(fit8stands_standpred)
for(i in 1:4){
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = 0,
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white') 
  
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = time[i,'warmup'],
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup']+time[i,'sample'],
       col =  '#87C55F',
       border = 'white') 
}

k <- 4
time <- get_elapsed_time(fit8stands_standpred_improved)
for(i in 1:4){
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = 0,
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white') 
  
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = time[i,'warmup'],
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup']+time[i,'sample'],
       col =  '#87C55F',
       border = 'white') 
}

k <- 5
time <- fit8stands_standpred_improveds_cmdstand$time()$chains
for(i in 1:4){
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = 0,
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white') 
  
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = time[i,'warmup'],
       xright = 0.5+i+k*5,
       ytop =  time[i,'warmup']+time[i,'sampling'],
       col =  '#87C55F',
       border = 'white') 
}
