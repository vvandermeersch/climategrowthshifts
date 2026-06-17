

round_superior <- function(x) {
  ceiling(x * 10) / 10
}

round_inferior <- function(x) {
  floor(x * 10) / 10
}

pdf('/home/victor/projects/climategrowthshifts/docs/misc/egu26/figures/climvar.pdf',
    height = 14 / 2.54, width = 35 / 2.54)
par(mar = c(4.5,1,1,1), cex.axis = 1.8, cex.lab = 2)
plot(x = NULL, y = NULL,
     xlim = c(-5, 7),
     ylim = c(0, 1.3),
     yaxt = 'n', xaxt = 'n',
     xlab = 'Average relative change in growth due to climate predictors (%)',
     ylab = '',
     main = '',
     frame.plot	= FALSE)

axis(1, at = seq(-5, 7, 1),
     labels = seq(-5, 7, 1),,
     las = 1)


for(i in 1:data$N_all_years){
  x <- gq_samples[[paste0('avg_exp_m1_clim[',i,']')]]
  param <- hist(x, breaks = seq(round_inferior(range(x)[1]), round_superior(range(x)[2]), 0.1), plot = FALSE)
  param_counts <- param$counts / max(param$counts) 
  rect(xleft = param$breaks[1:(length(param$breaks)-1)],
       ybottom = 0,
       xright = param$breaks[2:(length(param$breaks))],
       ytop =  param_counts ,
       col = 'grey90',
       border = 'white')
}

x <- gq_samples[[paste0('avg_exp_m1_clim_period[','30s',']')]]
param <- hist(x, breaks = seq(round_inferior(range(x)[1]), round_superior(range(x)[2]), 0.1), plot = FALSE)
param_counts <- param$counts / max(param$counts) 
rect(xleft = param$breaks[1:(length(param$breaks)-1)],
     ybottom = 0,
     xright = param$breaks[2:(length(param$breaks))],
     ytop =  param_counts ,
     col = '#cf6448',
     border = 'white')

x <- gq_samples[[paste0('avg_exp_m1_clim_period[','post2000',']')]]
param <- hist(x, breaks = seq(round_inferior(range(x)[1]), round_superior(range(x)[2]), 0.1), plot = FALSE)
param_counts <- param$counts / max(param$counts) 
rect(xleft = param$breaks[1:(length(param$breaks)-1)],
     ybottom = 0,
     xright = param$breaks[2:(length(param$breaks))],
     ytop =  param_counts ,
     col = '#cf6448',
     border = 'white')

x <- gq_samples[[paste0('avg_exp_m1_clim_period[','80s',']')]]
param <- hist(x, breaks = seq(round_inferior(range(x)[1]), round_superior(range(x)[2]), 0.1), plot = FALSE)
param_counts <- param$counts / max(param$counts) 
rect(xleft = param$breaks[1:(length(param$breaks)-1)],
     ybottom = 0,
     xright = param$breaks[2:(length(param$breaks))],
     ytop =  param_counts ,
     col = '#50938a',
     border = 'white')

x <- gq_samples[[paste0('avg_exp_m1_clim_period[','90s',']')]]
param <- hist(x, breaks = seq(round_inferior(range(x)[1]), round_superior(range(x)[2]), 0.1), plot = FALSE)
param_counts <- param$counts / max(param$counts) 
rect(xleft = param$breaks[1:(length(param$breaks)-1)],
     ybottom = 0,
     xright = param$breaks[2:(length(param$breaks))],
     ytop =  param_counts ,
     col = '#50938a',
     border = 'white')
dev.off()