inv_logit <- function(x){
  return(exp(x)/(1+exp(x)))
}


mu <- rnorm(1e5, 4.5, 0.5)
plot(density(inv_logit(mu)), col = 'darkblue', xlim = c(0.8,1))
tau <- abs(rnorm(1e5, 0, 1))
alpha <- rnorm(1e5, mu, tau)
y <- inv_logit(alpha)
lines(density(y), col = 'lightblue')
abline(v = quantile(inv_logit(mu), c(0.05, 0.95)), lty = 2, col = 'darkblue')
abline(v = quantile(y, c(0.05, 0.95)), lty = 2, col = 'lightblue')

par(mfrow = c(2,1))
phi_prior <- rbeta(1e6, 10, 500)
hist(phi_prior, breaks=seq(0, 1, 0.001),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab=expression(phi), yaxt='n', ylab="", xlim = c(0,0.25),
     mar = c(1,1,1,1))
quantile(phi_prior, c(0.05, 0.95))
abline(v = quantile(phi_prior, c(0.05, 0.95)), lty = 2, col = 'black')
omega_prior <- rbeta(1e6, 5, 5)
hist(omega_prior, breaks=seq(0, 1, 0.001),
     main="", col=util$c_light, border=util$c_light_highlight,
     xlab=expression(omega), yaxt='n', ylab="", xlim = c(0.,1),
     mar = c(1,1,1,1))
quantile(omega_prior, c(0.05, 0.95))
abline(v = quantile(omega_prior, c(0.05, 0.95)), lty = 2, col = 'black')
