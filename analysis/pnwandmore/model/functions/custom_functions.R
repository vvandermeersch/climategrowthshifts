plot_disc_pushforward_quantiles_2clades <- function(samples, names,
                                                    clades,
                                            baseline_values=NULL,
                                            baseline_col="black",
                                            residual=FALSE,
                                            xlab="", xticklabs=NULL,
                                            ylab=NULL, display_ylim=NULL, yticklabs=NULL,
                                            main="", line0 = FALSE) {
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
         ylim=display_ylim, ylab=ylab,  frame.plot=F)
    if(line0){abline(h = 0, lty = 2, col = "grey70")}
  } else {
    if (length(xticklabs) == N & is.null(yticklabs)) {
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab, xaxt="n",
           ylim=display_ylim, ylab=ylab,  frame.plot=F)
      if(line0){abline(h = 0, lty = 2, col = "grey70")}
      axis(1, at=1:N, labels=xticklabs)
    } else if(is.na(xticklabs)  & is.null(yticklabs)){
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab, xaxt="n",
           ylim=display_ylim, ylab=ylab,  frame.plot=F)
      if(line0){abline(h = 0, lty = 2, col = "grey70")}
    } else if(is.na(xticklabs)  & is.na(yticklabs) ){
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab, xaxt="n",
           ylim=display_ylim, ylab=ylab,  yaxt = "n",  frame.plot=F)
      if(line0){abline(h = 0, lty = 2, col = "grey70")}
    } else {
      warning(paste0('The list of x labels tick has the wrong',
                     ' dimension and baselines will not be plotted.'))
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab,
           ylim=display_ylim, ylab=ylab,  frame.plot=F)
      if(line0){abline(h = 0, lty = 2, col = "grey70")}
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



plot_expectand_pushforward_reverse <- function(expectand_vals, B,
         display_name="f",
         flim=NULL, ylim=NULL,
         col=c_dark, border="#DDDDDD",
         add=FALSE, main="", 
         baseline=NULL,
         baseline_col="black") {
  validate_array(expectand_vals, 'expectand_vals')
  
  # Automatically adjust histogram range to range of expectand values
  # if range is not already set as an input variable
  if (is.null(flim)) {
    min_f <- min(expectand_vals)
    max_f <- max(expectand_vals)
    delta <- (max_f - min_f) / B
    
    # Add bounding bins
    B <- B + 2
    min_f <- min_f - delta
    max_f <- max_f + delta
    flim <- c(min_f, max_f)
    
    bins <- seq(min_f, max_f, delta)
  } else {
    min_f <- flim[1]
    max_f <- flim[2]
    
    delta <- (max_f - min_f) / B
    bins <- seq(min_f, max_f, delta)
  }
  
  # Check value containment
  S <- dim(expectand_vals)[1] * dim(expectand_vals)[2]
  
  S_low <- sum(c(expectand_vals, recursive=TRUE) < min_f)
  if (S_low == 1)
    warning(sprintf('%i value (%.1f%%) fell below the histogram binning.',
                    S_low, 100 * S_low / S))
  else if (S_low > 1)
    warning(sprintf('%i values (%.1f%%) fell below the histogram binning.',
                    S_low, 100 * S_low / S))
  
  S_high <- sum(max_f < c(expectand_vals, recursive=TRUE))
  if (S_low == 1)
    warning(sprintf('%i value (%.1f%%) fell above the histogram binning.',
                    S_high, 100 * S_high / S))
  else if (S_low > 1)
    warning(sprintf('%i values (%.1f%%) fell above the histogram binning.',
                    S_high, 100 * S_high / S))
  
  # Compute bin heights
  mean_p <- rep(0, B)
  delta_p <- rep(0, B)
  
  for (b in 1:B) {
    # Estimate bin probabilities
    bin_indicator <- function(x) {
      ifelse(bins[b] <= x & x < bins[b + 1], 1, 0)
    }
    indicator_vals <- util$eval_uni_expectand_pushforward(expectand_vals,
                                                     bin_indicator)
    est <- util$ensemble_mcmc_est(indicator_vals)
    
    # Normalize bin probabilities by bin width to allow
    # for direct comparison to probability density functions
    width = bins[b + 1] - bins[b]
    mean_p[b] = est[1] / width
    delta_p[b] = est[2] / width
  }
  
  # Plot histogram
  idx <- rep(1:B, each=2)
  x <- sapply(1:length(idx), function(b) if(b %% 2 == 1) bins[idx[b]]
              else bins[idx[b] + 1])
  lower_inter <- sapply(idx, function (n)
    max(mean_p[n] - 2 * delta_p[n], 0))
  upper_inter <- sapply(idx, function (n)
    min(mean_p[n] + 2 * delta_p[n], 1 / width))
  
  if (add) {
    polygon(c(lower_inter, rev(upper_inter)), c(x, rev(x)),
            col=border, border=NA)
    lines(mean_p[idx], x,  col=col, lwd=2)
  } else {
    if (is.null(ylim)) {
      ylim=c(0, max(1.05 * upper_inter))
    }
    
    plot(1, type="n", main=main,
         xlim=ylim, xlab=display_name, 
         ylim=flim, ylab="", xaxt="n", frame.plot=F)
    # title(ylab="Estimated Bin\nProbabilities / Bin Width",
    #       mgp=c(1, 1, 0))
    
    polygon(c(lower_inter, rev(upper_inter)), c(x, rev(x)),
            col=border, border=NA)
    lines(mean_p[idx], x,  col=col, lwd=2)
  }
  
  # Plot baseline if applicable
  if (!is.null(baseline)) {
    abline(v=baseline, col="white", lty=1, lwd=4)
    abline(v=baseline, col=baseline_col, lty=1, lwd=2)
  }
}
