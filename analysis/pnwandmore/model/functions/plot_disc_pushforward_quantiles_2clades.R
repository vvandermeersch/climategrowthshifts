plot_disc_pushforward_quantiles_2clades <- function(samples, names,
                                                    clades,
                                            baseline_values=NULL,
                                            baseline_col="black",
                                            residual=FALSE,
                                            xlab="", xticklabs=NULL,
                                            ylab=NULL, display_ylim=NULL, yticklabs=NULL,
                                            main="") {
  # Check that baseline values are well-defined
  if (!is.null(baseline_values)) {
    if (length(baseline_values) != length(names)) {
      warning(paste0('The list of baseline values has the wrong',
                     ' dimension.  Baselines will not be plotted.'))
      baseline_values <- NULL
    }
  }
  
  # Check that names are in samples
  names <- check_expectand_names(names, samples)
  
  # Construct bins
  N <- length(names)
  bin_min <- 0.5
  bin_max <- N + 0.5
  bin_delta <- 1
  breaks <- seq(bin_min, bin_max, bin_delta)
  
  plot_config <- configure_bin_plotting(breaks)
  plot_idxs <- plot_config[[1]]
  plot_xs <- plot_config[[2]]
  
  # Construct marginal quantiles
  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  
  if(!is.null(baseline_values) & residual) {
    calc <- function(n) {
      util$ensemble_mcmc_quantile_est(samples[[names[n]]] -
                                        baseline_values[n],
                                      probs)
    }
    quantiles <- sapply(1:N, calc)
  } else {
    calc <- function(n) {
      util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
    }
    quantiles <- sapply(1:N, calc)
  }
  
  plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                          function(n) quantiles[1:9, n]))
  
  # Plot
  if (is.null(display_ylim)) {
    if (is.null(baseline_values)) {
      display_ylim <- c(min(quantiles[1,]),
                        max(quantiles[9,]))
    }
    else {
      if (residual) {
        display_ylim <- c(min(c(0, quantiles[1,])),
                          max(c(0, quantiles[9,])))
      } else {
        display_ylim <- c(min(min(quantiles[1,]),
                              min(baseline_values)),
                          max(max(quantiles[9,]),
                              max(baseline_values)))
      }
    }
    delta <- 0.05 * (display_ylim[2] - display_ylim[1])
    display_ylim[1] <- display_ylim[1] - delta
    display_ylim[2] <- display_ylim[2] + delta
  }
  
  if (is.null(ylab)) {
    if (is.null(baseline_values) | !residual)
      ylab <- "Marginal Quantiles"
    else
      ylab <- "Marginal Quantiles - Baselines"
  }
  
  if (is.null(xticklabs)) {
    plot(1, type="n", main=main,
         xlim=c(bin_min, bin_max), xlab=xlab,
         ylim=display_ylim, ylab=ylab)
  } else {
    if (length(xticklabs) == N & is.null(yticklabs)) {
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab, xaxt="n",
           ylim=display_ylim, ylab=ylab)
      axis(1, at=1:N, labels=xticklabs)
    } else if(is.na(xticklabs)  & is.na(yticklabs) ){
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab, xaxt="n",
           ylim=display_ylim, ylab=ylab,  yaxt = "n")
    } else {
      warning(paste0('The list of x labels tick has the wrong',
                     ' dimension and baselines will not be plotted.'))
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab,
           ylim=display_ylim, ylab=ylab)
    }
  }
  
  
  for (n in 1:N) {
    idx1 <- 2 * n - 1
    idx2 <- 2 * n
    
    polygon(c(plot_xs[idx1:idx2], rev(plot_xs[idx1:idx2])),
            c(plot_quantiles[1,idx1:idx2], rev(plot_quantiles[9,idx1:idx2])),
            col = ifelse(clades[n] == 1, '#9393c7', '#93c793'), border = NA)
    polygon(c(plot_xs[idx1:idx2], rev(plot_xs[idx1:idx2])),
            c(plot_quantiles[2,idx1:idx2], rev(plot_quantiles[8,idx1:idx2])),
            col = ifelse(clades[n] == 1, '#7d7dbb', '#7dbb7d'), border = NA)
    polygon(c(plot_xs[idx1:idx2], rev(plot_xs[idx1:idx2])),
            c(plot_quantiles[3,idx1:idx2], rev(plot_quantiles[7,idx1:idx2])),
            col = ifelse(clades[n] == 1, '#6767b0', '#67b067'), border = NA)
    polygon(c(plot_xs[idx1:idx2], rev(plot_xs[idx1:idx2])),
            c(plot_quantiles[4,idx1:idx2], rev(plot_quantiles[6,idx1:idx2])),
            col = ifelse(clades[n] == 1, '#5252a5', '#52a552'), border = NA)
    
    lines(plot_xs[idx1:idx2], plot_quantiles[5, idx1:idx2],
          col=ifelse(clades[n] == 1, '#27278f', '#278f27'), lwd=2)
  }
  
  if (!is.null(baseline_values)) {
    if (residual) {
      abline(h=0, col="#DDDDDD", lwd=2, lty=3)
    } else {
      for (n in 1:N) {
        idx1 <- 2 * n - 1
        idx2 <- 2 * n
        lines(plot_xs[idx1:idx2], rep(baseline_values[n], 2),
              col="white", lwd=4)
        lines(plot_xs[idx1:idx2], rep(baseline_values[n], 2),
              col=baseline_col, lwd=2)
      }
    }
  }
}
