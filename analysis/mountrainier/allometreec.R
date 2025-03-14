
rm(list=ls())
wd <- "~/projects/climategrowthshifts/analysis/mountrainier"

library(paletteer)

gaussfunc = function(x, mu, sigma){
  return(exp(-((x - mu)^2)/(2*(sigma^2))))
}

dm <- function(t,m,mmax,mu){
  return(gaussfunc(t, mu, 60)* (m/mmax))
}

allometree <- function(r, a, b){
  return(exp(b*(log(r*2)) + a))
}

par(mfrow = c(1,2))
t <- seq(1,100,1)
plot(dm(seq(1,100,1), m = 100, mmax = 100, mu = 50) ~ t, cex = 0.5, , xlab = 'year', ylab = 'dm/dt')
r <- seq(2, 50,1)
plot(allometree(r, a = 0.2, b = 0.55) ~ c(r*2), cex = 0.5, , xlab = 'dbh (cm)', ylab = 'height (m)')


rho <- 700*1e-6 # kg.cm-3

forest <- data.frame()
for(i in 1:50){
  
  r0 <- runif(1, 5, 10)
  mmax <- runif(1, 80, 120)
  mu <- runif(1, 40, 60)
  
  a = runif(1, 0.1, 0.3)
  b = runif(1, 0.6, 0.7)
  
  t <- 1
  dr <- c()
  r <- c(r0)
  ba <- c(pi*r[t]^2)
  h <- allometree(r[t], a, b)*100 
  m <- c(rho*ba[t]*(h[t]))
  
  for(t in 2:200){
    
    dr[t-1] <- 1/sqrt(pi*rho*h[t-1]*m[t-1])*dm(t-1,m[t-1], mmax, mu)
    r[t] <- r[t-1] + dr[t-1]
    
    ba[t] <- pi*r[t]^2
    h[t] <-  allometree(r[t], a, b)*100 # allometric equation
    m[t] <- rho*ba[t]*(h[t])
    
    
  }
  
  # plot(r*2 ~ c(1:200), cex = 0.5, xlab = 'time', ylab = 'diameter')
  # plot(dr ~ c(1:199), cex = 0.5, xlab = 'time', ylab = 'ringwidth')
  # plot(dr ~ m[1:199], cex = 0.5)
  
  forest <- rbind(forest, data.frame(tree = i, time = 1:200, mass = m, dbh = r*2, height = h, rw = c(NA, dr)))
  
}


colors = paletteer_c("viridis::inferno", n=nColor)
rank <- as.factor( as.numeric( cut(forest$tree, 20)))


pdf(file.path(wd, 'figures', 'allometric_simulations.pdf'), width = 8, height = 7)
par(mfrow = c(2,2))
plot(forest$dbh ~ forest$time, cex = 0.1, 
     col = colors[ rank ],
     xlab = 'year', ylab = 'diameter (cm)')

plot(forest$rw*10 ~ forest$time, cex = 0.1, 
     col = colors[ rank ],
     xlab = 'year', ylab = 'ringwidth (mm)')

plot(forest$rw*10 ~ forest$mass, cex = 0.1, 
     col = colors[ rank ],
     xlab = 'mass (kg)', ylab = 'ringwidth (mm)')

plot(forest$h/100 ~ forest$dbh, cex = 0.05, 
     col = colors[ rank ],
     xlab = 'dbh (cm)', ylab = 'height (m)')
dev.off()
