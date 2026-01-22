
logit <- function(p){
  return(log(p/(1-p)))
}

invlogit <- function(x){
  return(1/(1+exp(-x)))
}

phi_sck0 <- rbeta(1,2 ,20)
pre0 <- 5
pre <- seq(0,20,0.1)
beta_shock_pre <- -0.1

plot(invlogit(logit(phi_sck0) + beta_shock_pre * (pre-pre0)) ~ pre, type = 'l')
abline(v = pre0, lty = 2)
abline(v = 2*pre0, lty = 2)

phi_sck0 <- rbeta(1e6,2 ,20)
hist(phi_sck0)

# with beta = -0.15, doubling precipitation divide the probability by 2:
invlogit(logit(0.1) - 0.15 * (5 - 5))
invlogit(logit(0.1) - 0.15 * (10 - 5))
# so a prior like normal(0, 0.15 / 2.57)?