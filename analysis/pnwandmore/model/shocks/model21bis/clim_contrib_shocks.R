
s <- 1
years <- seq(1+data$N_all_years*(s-1), data$N_all_years*s, 1)

pre2000 <- years[which(data$all_years[years] < 2000)]
post2000 <- years[which(data$all_years[years] >= 2000)]

delta_vpd <- mean(data$vpd_obs[post2000]) - mean(data$vpd_obs[pre2000])
delta_pre <- mean(data$pre_obs[post2000]) - mean(data$pre_obs[pre2000])
delta_gdd <- mean(data$gdd_obs[post2000]) - mean(data$gdd_obs[pre2000])

contrib_vpd <- base_samples[['beta_phi_vpd']] * delta_vpd
contrib_pre <- base_samples[['beta_phi_pre']] * delta_pre
contrib_gdd <- base_samples[['beta_phi_gdd']] * delta_gdd

total <- contrib_vpd + contrib_pre + contrib_gdd 
share_vpd <- 100 * contrib_vpd / total
share_pre <- 100 * contrib_pre / total
share_gdd <- 100 * contrib_gdd / total

util$ensemble_mcmc_quantile_est(share_vpd, c(0.025, 0.5, 0.975))
util$ensemble_mcmc_quantile_est(share_pre, c(0.025, 0.5, 0.975))
util$ensemble_mcmc_quantile_est(share_gdd, c(0.025, 0.5, 0.975))


par(mfrow = c(1,1), mar = c(4,4,1,1))
util$plot_expectand_pushforward(share_pre, flim = c(0, 100), 20, col = '#4F6D7A', ylim = c(0, 0.15),
                                display_name = ' Contribution (%)')
util$plot_expectand_pushforward(share_gdd, flim = c(0, 100), 20, col = '#7A9E7E', add = T, ylim = c(0, 0.1))
util$plot_expectand_pushforward(share_vpd, flim = c(0, 100), 20, col = '#C1665A', ylim = c(0, 0.15), add = T)
legend(
  "topright", inset = c(0.1, 0.1),
  legend = c('Winter precipitation', 'Summer VPD', 'Year-round GDD'),
  col = c('#4F6D7A', "#C1665A", '#7A9E7E'),
  pch = 15, cex = 1, bty = 'n', xpd = NA, pt.cex = 1.5
)


par(mfrow = c(1,1), mar = c(4,4,1,1))
util$plot_pairs_by_chain(base_samples[['beta_phi_vpd']], 'beta_phi_vpd',
                         base_samples[['beta_phi_gdd']], 'beta_phi_gdd')


util$plot_pairs_by_chain(share_gdd, 'Contribution GDD (%)',
                         share_vpd, 'Contribution VPD (%)')


