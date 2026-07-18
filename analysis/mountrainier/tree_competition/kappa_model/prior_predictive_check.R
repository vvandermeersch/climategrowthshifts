prior_y0 <- function(n=100,
                     mu_log_y0_mu,
                     tau_log_y0_mu,
                     tau_log_y0_tau){
  log_y0_mu <- rnorm(n,mu_log_y0_mu,tau_log_y0_mu)
  log_y0_tau <- rnorm(n,0,tau_log_y0_tau)
  y0_mu <- numeric(0)
  y0_var <- numeric(0)
  
  for (i in 1:n) {
    y0_raw <- rnorm(n)
    y0 <- exp(log_y0_mu[i]+log_y0_tau[i]*y0_raw)
    y0_mu <- c(y0_mu,mean(y0))
    y0_var <- c(y0_var,var(y0))
  }
  
  total_mu <- mean(y0_mu)
  total_var <- mean(y0_var)+var(y0_mu)
  return(
    c(
      mean=total_mu,
      var = total_var,
      sd = sqrt(total_var)
    )
  )
}

prior_r <- function(n=100,
                    mu_log_odds_mu,
                    tau_log_odds_mu,
                    tau_log_odds_tau) {
  log_odds_mu <- rnorm(n,mu_log_odds_mu,tau_log_odds_mu)
  log_odds_tau <- abs(rnorm(n,0,tau_log_odds_tau))
  r_mean <- numeric(n)
  r_var <- numeric(n)
  
  for (i in 1:n) {
    log_odds_raw <- rnorm(n)
    log_odds <- log_odds_mu[i]+log_odds_raw*log_odds_tau[i]
    r <- plogis(log_odds)
    
    r_mean[i] <- mean(r)
    r_var[i] <- var(r)
  }
  
  return(c(mean = mean(r_mean),
           var = mean(r_var)+var(r_mean)))
}
