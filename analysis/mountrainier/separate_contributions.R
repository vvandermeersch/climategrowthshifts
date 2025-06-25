# Posterior Quantification
fit <- stan(file=file.path(wd, 'stan/model7_wpredictors_sepcontrib.stan'),
            data=data, seed=5838299, cores = 4,
            warmup=1000, iter=2024, refresh=10)

samples <- util$extract_expectand_vals(fit)

t <- 1

nmin <- data$tree_start_idxs[t]
nmax <- data$tree_end_idxs[t]

pred_names <- sapply(nmin:nmax, function(n) paste0('fsum[', n, ']'))
N <- length(pred_names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[pred_names[n]]], probs)
}
fsum <- sapply(1:N, calc)

pred_names <- sapply(nmin:nmax, function(n) paste0('f[', n, ']'))
N <- length(pred_names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[pred_names[n]]], probs)
}
f <- sapply(1:N, calc)

pred_names <- sapply(nmin:nmax, function(n) paste0('mu[', n, ']'))
N <- length(pred_names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[pred_names[n]]], probs)
}
mu <- sapply(1:N, calc)

pred_names <- sapply(nmin:nmax, function(n) paste0('log_rw[', n, ']'))
N <- length(pred_names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[pred_names[n]]], probs)
}
log_rw <- sapply(1:N, calc)

x <- years[nmin:nmax]

ggplot() +
  geom_ribbon(aes(x = x, ymin = fsum[1,], ymax = fsum[9,]), fill =  '#bebedd', alpha = 0.5) +
  #geom_ribbon(aes(x = x, ymin = fsum[2,], ymax = fsum[8,]), fill =  '#a8a8d2', alpha = 0.5) +
  #geom_ribbon(aes(x = x, ymin = fsum[3,], ymax = fsum[7,]), fill =  '#9393c7', alpha = 0.5) +
  #geom_ribbon(aes(x = x, ymin = fsum[4,], ymax = fsum[6,]), fill =  '#7d7dbb', alpha = 0.5) +
  geom_ribbon(aes(x = x, ymin = f[1,], ymax = f[9,]), fill =  '#a4cfa4', alpha = 0.5) +
  #geom_ribbon(aes(x = x, ymin = f[2,], ymax = f[8,]), fill =  '#97c897', alpha = 0.5) +
  #geom_ribbon(aes(x = x, ymin = f[3,], ymax = f[7,]), fill =  '#8ac18a', alpha = 0.5) +
  #geom_ribbon(aes(x = x, ymin = f[4,], ymax = f[6,]), fill =  '#7dbb7d', alpha = 0.5) +
  
  geom_line(aes(x = x, y = log_rw[5,]), color =  'darkorange', size = 0.6) +
  # geom_line(aes(x = x, y = mu[5,]), color =  'darkred', size = 0.6) +
  geom_line(aes(x = x, y = f[5,]), color =  '#7dbb7d', linetype = 'dashed', size = 1) +
  geom_line(aes(x = x, y = fsum[5,]), color =  '#7d7dbb', linetype = 'dotted', size = 1.5) +
  
  
  theme_bw() +
  annotate('label', x = 1980, y = -0.18, color = '#7d7dbb', label = 'mu1+mu2+mu3', hjust = 0, vjust = 0, size = 5) +
  annotate('label', x = 1983, y = 0.32, color = '#7dbb7d', label = 'f', hjust = 0, vjust = 0, size = 7) +
  # annotate('label', x = 2004, y = 0.34, color = 'darkred', label = 'mu', hjust = 0, vjust = 0, size = 5) +
  annotate('label', x = 2000, y = 0.37, color = 'darkorange', label = 'log_rw', hjust = 0, vjust = 0, size = 5) +
  labs(y = '', x = 'Year')


pred_names <- sapply(nmin:nmax, function(n) paste0('mu2[', n, ']'))
N <- length(pred_names)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[pred_names[n]]], probs)
}
mu1 <- sapply(1:N, calc)

ggplot() +
  geom_ribbon(aes(x = x, ymin = mu1[1,], ymax = mu1[9,]), fill =  '#bebedd', alpha = 0.5) +
  geom_ribbon(aes(x = x, ymin = mu1[2,], ymax = mu1[8,]), fill =  '#a8a8d2', alpha = 0.5) +
  geom_ribbon(aes(x = x, ymin = mu1[3,], ymax = mu1[7,]), fill =  '#9393c7', alpha = 0.5) +
  geom_ribbon(aes(x = x, ymin = mu1[4,], ymax = mu1[6,]), fill =  '#7d7dbb', alpha = 0.5) +
  theme_bw() +
  labs(y = '', x = 'Year')
