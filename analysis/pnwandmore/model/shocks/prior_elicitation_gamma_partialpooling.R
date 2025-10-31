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
