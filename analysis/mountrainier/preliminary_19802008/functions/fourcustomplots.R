
fourcustomplots <- function(idx, samples, data, x, xbreaks = seq(1985, 2005,5), xannot = 1981, xname = ""){
  
  x <- x[idxs]
  y <- data$log_rw_obs[idxs]
  
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_ltonly_pred[', n, ']'))
  N <- length(rw_names)
  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  calc <- function(n) {
    util$ensemble_mcmc_quantile_est(samples[[rw_names[n]]], probs)
  }
  plot_quantiles_ltonly <- sapply(1:N, calc)
  
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_woclim_pred[', n, ']'))
  N <- length(rw_names)
  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  calc <- function(n) {
    util$ensemble_mcmc_quantile_est(samples[[rw_names[n]]], probs)
  }
  plot_quantiles_woclim <- sapply(1:N, calc)
  
  rw_names <- sapply(idxs,
                     function(n) paste0('log_rw_pred[', n, ']'))
  N <- length(rw_names)
  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  calc <- function(n) {
    util$ensemble_mcmc_quantile_est(samples[[rw_names[n]]], probs)
  }
  plot_quantiles <- sapply(1:N, calc)
  
  
  ltonly <- ggplot() +
    geom_ribbon(aes(x = x, ymin = plot_quantiles_ltonly[1,], ymax = plot_quantiles_ltonly[9,]), fill =  '#bebedd') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles_ltonly[2,], ymax = plot_quantiles_ltonly[8,]), fill =  '#a8a8d2') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles_ltonly[3,], ymax = plot_quantiles_ltonly[7,]), fill =  '#9393c7') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles_ltonly[4,], ymax = plot_quantiles_ltonly[6,]), fill =  '#7d7dbb') +
    geom_point(aes(x = x, y = y), size = 2, color = 'white') +
    geom_point(aes(x = x, y = y), size = 1, color = 'black') +
    theme_bw() +
    scale_x_continuous(expand = c(0,0), breaks = xbreaks) +
    scale_y_continuous(expand = c(0,0), limits = c(-2,2)) + 
    labs(x = xname, y = "log(ring width)") +
    annotate('text', x = xannot, y = Inf, color = '#7d7dbb', label = 'Long-term scale', hjust = 0, vjust = 2)
  
  woclim <- ggplot() +
    geom_ribbon(aes(x = x, ymin = plot_quantiles_woclim[1,], ymax = plot_quantiles_woclim[9,]), fill =  '#a4cfa4') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles_woclim[2,], ymax = plot_quantiles_woclim[8,]), fill =  '#97c897') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles_woclim[3,], ymax = plot_quantiles_woclim[7,]), fill =  '#8ac18a') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles_woclim[4,], ymax = plot_quantiles_woclim[6,]), fill =  '#7dbb7d') +
    geom_point(aes(x = x, y = y), size = 2, color = 'white') +
    geom_point(aes(x = x, y = y), size = 1, color = 'black') +
    theme_bw() +
    scale_x_continuous(expand = c(0,0), breaks = xbreaks) +
    scale_y_continuous(expand = c(0,0), limits = c(-2,2)) +
    labs(x = xname, y = "log(ring width)") +
    annotate('text', x = xannot, y = Inf, color = '#7dbb7d', label = 'Long- and short-term scales', hjust = 0, vjust = 2)
  
  full <- ggplot() +
    geom_ribbon(aes(x = x, ymin = plot_quantiles[1,], ymax = plot_quantiles[9,]), fill =  '#cfa4a4') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles[2,], ymax = plot_quantiles[8,]), fill =  '#c89797') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles[3,], ymax = plot_quantiles[7,]), fill =  '#c18a8a') +
    geom_ribbon(aes(x = x, ymin = plot_quantiles[4,], ymax = plot_quantiles[6,]), fill =  '#bb7d7d') +
    geom_point(aes(x = x, y = y), size = 2, color = 'white') +
    geom_point(aes(x = x, y = y), size = 1, color = 'black') +
    theme_bw() +
    scale_x_continuous(expand = c(0,0), breaks = xbreaks) +
    scale_y_continuous(expand = c(0,0), limits = c(-2,2)) +
    labs(x = xname, y = "log(ring width)") +
    annotate('text', x = xannot, y = Inf, color = '#bb7d7d', label = 'Full model', hjust = 0, vjust = 2)
  
  all <- ggplot() +
    geom_line(aes(x = x, y = plot_quantiles_ltonly[5,]), color =  '#7d7dbb', linewidth = 1) + 
    geom_line(aes(x = x, y = plot_quantiles_woclim[5,]), color =  '#7dbb7d', linewidth = 1) + 
    geom_line(aes(x = x, y = plot_quantiles[5,]), color =  '#bb7d7d', linewidth = 1) + 
    geom_point(aes(x = x, y = y), size = 2, color = 'white') +
    geom_point(aes(x = x, y = y), size = 1, color = 'black') +
    theme_bw() +
    scale_x_continuous(expand = c(0,0), breaks = xbreaks) +
    scale_y_continuous(expand = c(0.15,0.15)) +
    labs(x = xname, y = "log(ring width)") +
    annotate('text', x = xannot, y = Inf, color = 'black', label = 'Zoom-in', hjust = 0, vjust = 2)
  
  
  plots <- ltonly + woclim + full + all + plot_layout(axis_titles = "collect", nrow = 2)
  
  return(plots)
}


