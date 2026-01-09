wd <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore'
library(rstan)

fit8stands <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO.rds'))
fit16stands <- readRDS(file.path(wd, 'output/model/scaling', 'fit_16standsPIPO.rds'))
fit8stands_ov <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO_oldversion.rds'))
fit8stands_standp <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO_standpred.rds'))

plot(x = NULL,
     y = NULL,
     xlim = c(0, 15),
     ylim = c(0, 36000),
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = 'Runtime (hours)',
     bty = "n")
par(mgp = c(3, 1.5, 0), mfrow = c(1,1))

axis(1, at = c(2.5,7.5, 12.5), labels = c('8 stands\n(initial version)', '8 stands\n(4 cores)', '8 stands\n(stand, 4 cores)'))
axis(2, at = seq(0, 36000, 3600), labels = sapply(seq(0, 10, 1), function(i) ifelse(i%%2==0,i,'')))

time <- get_elapsed_time(fit8stands_ov)
for(i in 1:4){
  rect(xleft = 0.5+(i-1),
       ybottom = 0,
       xright = 0.5+i,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white') 
  
  rect(xleft = 0.5+(i-1),
       ybottom = time[i,'warmup'],
       xright = 0.5+i,
       ytop =  time[i,'warmup']+time[i,'sample'],
       col =  '#87C55F',
       border = 'white') 
}

time <- get_elapsed_time(fit8stands)
for(i in 1:4){
  rect(xleft = 5.5+(i-1),
       ybottom = 0,
       xright = 5.5+i,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white')
  rect(xleft = 5.5+(i-1),
       ybottom = time[i,'warmup'],
       xright = 5.5+i,
       ytop =  time[i,'warmup']+time[i,'sample'],
       col =  '#87C55F',
       border = 'white') 
}

time <- get_elapsed_time(fit8stands_standp)
for(i in 1:4){
  rect(xleft = 10.5+(i-1),
       ybottom = 0,
       xright = 10.5+i,
       ytop =  time[i,'warmup'],
       col =  '#F6CF71',
       border = 'white')
  rect(xleft = 10.5+(i-1),
       ybottom = time[i,'warmup'],
       xright = 10.5+i,
       ytop =  time[i,'warmup']+time[i,'sample'],
       col =  '#87C55F',
       border = 'white') 
}
