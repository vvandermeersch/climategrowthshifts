
plot(x = NULL,
     y = NULL,
     xlim = c(0, 15),
     ylim = c(0, 0.1),
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = 'Stepsize',
     bty = "n")
par(mgp = c(3, 1.5, 0))

axis(1, at = c(2.5,7.5, 12.5), labels = c('8 stands\n(initial version)', '8 stands\n(4 cores)', '16 stands\n(4 cores)'))
axis(2, at = seq(0, 0.1, 0.02), labels = seq(0, 0.1, 0.02))

params <- get_sampler_params(fit8stands_ov, inc_warmup = FALSE)
for(i in 1:4){
  start <- (1+(i-1)*1000)
  end <- start+999
  rect(xleft = 0.5+(i-1),
       ybottom = 0,
       xright = 0.5+i,
       ytop =  unique(params[[i]][,2]),
       col =  '#F89C74',
       border = 'white') 
  
}

params <- get_sampler_params(fit8stands, inc_warmup = FALSE)
for(i in 1:4){
  start <- (1+(i-1)*1000)
  end <- start+999
  rect(xleft = 5.5+(i-1),
       ybottom = 0,
       xright = 5.5+i,
       ytop =  unique(params[[i]][,2]),
       col =  '#F89C74',
       border = 'white') 
  
}

params <- get_sampler_params(fit16stands, inc_warmup = FALSE)
for(i in 1:4){
  start <- (1+(i-1)*1000)
  end <- start+999
  rect(xleft = 10.5+(i-1),
       ybottom = 0,
       xright = 10.5+i,
       ytop =  unique(params[[i]][,2]),
       col =  '#F89C74',
       border = 'white') 
  
}