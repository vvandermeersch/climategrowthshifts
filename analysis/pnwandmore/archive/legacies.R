





n <- 1
R <- matrix(c(1, 0.7, 0.25,
              0.7, 1, 0.7,
              0.25, 0.7, 1), 
            nrow = 3, ncol = 3, byrow = TRUE)

mu <- c(X = 0, Y = 0, Z = 0)



MASS::mvrnorm(n, mu = mu, Sigma = R)

exp_quad <- function(dx, gamma, rho){
  return(gamma^2 * exp(-1/2*(dx/rho)^2))
}

plotdf <- data.frame(dist = seq(0, 4, 1))


beta0 <- -0.09

gamma = 1
rho <- 0.7
plotdf$cov <- sapply(plotdf$dist, function(i) exp_quad(i, gamma = 1, rho = rho))
plotdf$beta <- beta0*plotdf$cov 

par(mfrow = c(2,1))
plot.new()
plot.window(xlim = c(0,max(plotdf$dist)), ylim = c(0,gamma))
grid()
axis(1, las = 1, cex.axis = 0.8, tck=-0.02, at = c(0,rho, 5), labels=c("0", expression(rho), ""))
axis(2, las = 1, cex.axis = 0.8, tck=-0.02, at = c(0,gamma^2/exp(1/2),gamma^2), labels=c("0", expression(frac(gamma^2,exp(1/2))), expression(gamma^2)))
title(main = "Covariogram", line=1, cex.main = 1)
title(xlab = expression(paste(Delta,year)),  line=3, cex.main = 1)
lines(plotdf$dist, plotdf$cov, col="white", lwd=4)
lines(plotdf$dist, plotdf$cov, col="#4F5F28FF", lwd=2)

plot.new()
plot.window(xlim = c(0,max(plotdf$dist)), ylim = c(beta0,0))
grid()
axis(1, las = 1, cex.axis = 0.8, tck=-0.02, at = c(0,rho, 5), labels=c("0", expression(rho), ""))
axis(2, las = 1, cex.axis = 0.8, tck=-0.02, at = c(0,beta0*gamma^2/exp(1/2),beta0), labels=c("0", expression(frac(beta[0]%*%gamma^2,exp(1/2))), expression(beta[0])))
title(main = "Drought effect", line=1, cex.main = 1)
title(xlab = expression(paste(Delta,year)),  line=3, cex.main = 1)
lines(plotdf$dist, plotdf$beta, col="white", lwd=4)
lines(plotdf$dist, plotdf$beta, col="#AE2565FF", lwd=2)


