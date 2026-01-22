

fit_8stands <- readRDS(file.path(wd, 'output/model/scaling', 'fit_8standsPIPO.rds'))
fit_16stands <- readRDS(file.path(wd, 'output/model/scaling', 'fit_16standsPIPO.rds'))
fit_32stands <- readRDS(file.path(wd, 'output/model/scaling', 'fit_32standsPIPO_cmdstan.rds'))


plot(x = NULL,
     y = NULL,
     xlim = c(0, 15),
     ylim = c(0, 3600*24),
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = 'Runtime (hours)',
     bty = "n")
par(mgp = c(3, 1.5, 0), mfrow = c(1,1))

axis(1, at = seq(2.5,12.5,5), labels = c('8 stands', '16 stands', '32 stands'))
axis(2, at = seq(0, 3600*24, 3600*2), labels = sapply(seq(0, 24, 2), function(i) ifelse(i%%4==0,i,'')))

k <- 0
time <- get_elapsed_time(fit_8stands)
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
timeprev <- time

k <- 1
time <- get_elapsed_time(fit_16stands)
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
for(i in 1:4){
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = 0,
       xright = 0.5+i+k*5,
       ytop =  (timeprev[i,'warmup']+timeprev[i,'sample'])*2,
       col =  NA,
       border = 'black') 
}
timeprev <- time


k <- 2
time <- fit_32stands$time()$chains
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
for(i in 1:4){
  rect(xleft = 0.5+(i-1)+k*5,
       ybottom = 0,
       xright = 0.5+i+k*5,
       ytop =  (timeprev[i,'warmup']+timeprev[i,'sample'])*2,
       col =  NA,
       border = 'black') 
}

