
gaussfunc = function(x, mu, sigma){
  return(exp(-((x - mu)^2)/(2*(sigma^2))))
}

dm <- function(t,m,mmax,mu){
  
  return(gaussfunc(t, mu, 60)* (m/mmax))
  
}

allometree <- function(r){
  return(exp(0.5*(log(r*2)) + 0.2))
}

t <- seq(1,100,1)
plot(dm(seq(1,100,1), m = 100, mmax = 100, mu = 50) ~ t, cex = 0.5, , xlab = 'time', ylab = 'dm/dt, constant mass')

r <- seq(5,200,1)
plot(allometree(r) ~ c(r*2), cex = 0.5, , xlab = 'dbh', ylab = 'height')


rho <- 700*1e-6 # kg.cm-3

forest <- data.frame()
for(i in 1:100){
  
  r0 <- runif(1, 5, 10)
  mmax <- runif(1, 50, 200)
  mu <- runif(1, 40, 60)
  
  t <- 1
  dr <- c()
  r <- c(r0)
  ba <- c(pi*r[t]^2)
  h <- allometree(r[t])*100 
  m <- c(rho*ba[t]*(h[t]))
  
  for(t in 2:200){
    
    dr[t-1] <- 1/sqrt(pi*rho*h[t-1]*m[t-1])*dm(t-1,m[t-1], mmax, mu)
    r[t] <- r[t-1] + dr[t-1]
    
    ba[t] <- pi*r[t]^2
    h[t] <-  allometree(r[t])*100 # allometric equation
    m[t] <- rho*ba[t]*(h[t])
    
    
  }
  
  plot(r*2 ~ c(1:200), cex = 0.5, xlab = 'time', ylab = 'diameter')
  plot(dr ~ c(1:199), cex = 0.5, xlab = 'time', ylab = 'ringwidth')
  plot(dr ~ m[1:199], cex = 0.5)
  
  forest <- rbind(forest, data.frame(tree = i, time = 1:200, mass = m, dbh = r*2, height = h, rw = c(NA, dr)))
  
}

library(paletteer)
colors = paletteer_c("viridis::inferno", n=nColor)
rank <- as.factor( as.numeric( cut(forest$tree, 20)))

plot(forest$dbh ~ forest$time, cex = 0.1, 
     col = colors[ rank ],
     xlab = 'time', ylab = 'diameter')

plot(forest$rw ~ forest$time, cex = 0.1, 
     col = colors[ rank ],
     xlab = 'time', ylab = 'ringwidth')

plot(forest$rw ~ forest$mass, cex = 0.1, 
     col = colors[ rank ],
     xlab = 'mass', ylab = 'ringwidth')

plot(forest$h/100 ~ forest$dbh, cex = 0.1, 
     col = colors[ rank ],
     xlab = 'dbh', ylab = 'height')
