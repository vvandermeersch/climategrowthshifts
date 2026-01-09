




plot(x = NULL,
     y = NULL,
     xlim = c(0, 15),
     ylim = c(0, 300),
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = 'Number of leapfrog evaluations',
     bty = "n")
par(mgp = c(3, 1.5, 0))

axis(1, at = c(2.5,7.5, 12.5), labels = c('8 stands\n(initial version)', '8 stands\n(4 cores)', '16 stands\n(4 cores)'))
axis(2, at = seq(0, 300, 100), labels = seq(0, 300, 100))

lpfg <- rstan::get_num_leapfrog_per_iteration(fit8stands_ov)
for(i in 1:4){
  start <- (1+(i-1)*1000)
  end <- start+999
  rect(xleft = 0.5+(i-1),
       ybottom = 0,
       xright = 0.5+i,
       ytop =  mean(lpfg[start:end]),
       col =  '#F89C74',
       border = 'white') 

}

lpfg <- rstan::get_num_leapfrog_per_iteration(fit8stands)
for(i in 1:4){
  start <- (1+(i-1)*1000)
  end <- start+999
  rect(xleft = 5.5+(i-1),
       ybottom = 0,
       xright = 5.5+i,
       ytop =  mean(lpfg[start:end]),
       col =  '#F89C74',
       border = 'white') 
  
}

lpfg <- rstan::get_num_leapfrog_per_iteration(fit16stands)
for(i in 1:4){
  start <- (1+(i-1)*1000)
  end <- start+999
  rect(xleft = 10.5+(i-1),
       ybottom = 0,
       xright = 10.5+i,
       ytop =  mean(lpfg[start:end]),
       col =  '#F89C74',
       border = 'white') 
  
}



