rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnwandmore"
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('functions/custom_functions.R', local=util)
setwd(wd)
library(ggplot2)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_24sept2025.rds'))
samples <- readRDS(file.path(wd, 'output/model', 'samples_24sept2025_partialpooling2clades_gdd_amjjas_1interaction.rds'))

range(data$vpd_obs) # hPa
range(data$sm_obs) # % (m3.m-3)

vpd_sm <- unique(data.frame(vpd = data$vpd_obs, sm = data$sm_obs))
sm_breaks <- seq(9,40,2)
vpd_sm$sm_bin <- cut(vpd_sm$sm, breaks = sm_breaks)
midpoints <- (sm_breaks[-length(sm_breaks)] + sm_breaks[-1]) / 2
vpd_sm$sm_bin_mid <- midpoints[as.numeric(vpd_sm$sm_bin)]

count_df <- aggregate(vpd ~ sm_bin_mid, vpd_sm, FUN = length)
names(count_df)[2] <- 'count_points'
vpd_sm <- merge(vpd_sm, count_df)

ggplot(data = vpd_sm) +
  geom_boxplot(aes(x = sm_bin_mid, y = vpd, group = sm_bin_mid, col = count_points), outliers = FALSE) +
  scale_color_gradient(low = 'grey90', high = 'black') +
  geom_boxplot(aes(x = sm, y = 32), outliers = FALSE, width = 2) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Soil moisture (m3.m-3)', y = 'VPD (hPa)')

# baselines defined in the model
sm0 <- 25
vpd0 <- 8

alpha <- util$ensemble_mcmc_quantile_est(samples[['alpha']], probs = 0.5) # 0.065

beta_sm_angio <- util$ensemble_mcmc_quantile_est(samples[['mu_sm[2]']], probs = 0.5) # 0.065
beta_sm_gymno <- util$ensemble_mcmc_quantile_est(samples[['mu_sm[1]']], probs = 0.5) # 0.022

beta_vpd_angio <- util$ensemble_mcmc_quantile_est(samples[['mu_vpd[2]']], probs = 0.5) # -0.034
beta_vpd_gymno <- util$ensemble_mcmc_quantile_est(samples[['mu_vpd[1]']], probs = 0.5) # -0.071


low_sm <- as.numeric(quantile(data$sm_obs, 0.05)) # 19.10%
vpd <- c(seq(range(data$vpd_obs)[1], range(data$vpd_obs)[2], 0.1), low_vpd, high_vpd)
logrw <- alpha + beta_sm_angio*(low_sm-sm0) + beta_vpd_angio*(vpd-vpd0)
angio_lowsm <- data.frame(logrw, rw = exp(logrw), vpd, sm = low_sm, clade = "Angiosperms")
logrw <- alpha + beta_sm_gymno*(low_sm-sm0) + beta_vpd_gymno*(vpd-vpd0)
gymno_lowsm <- data.frame(logrw, rw = exp(logrw), vpd, sm = low_sm, clade = "Gymnosperms")

high_sm <- as.numeric(quantile(data$sm_obs, 0.95)) # 32.17%
vpd <- c(seq(range(data$vpd_obs)[1], range(data$vpd_obs)[2], 0.1), low_vpd, high_vpd)
logrw <- alpha + beta_sm_angio*(high_sm-sm0) + beta_vpd_angio*(vpd-vpd0)
angio_highsm <- data.frame(logrw, rw = exp(logrw), vpd, sm = high_sm, clade = "Angiosperms")
logrw <- alpha + beta_sm_gymno*(high_sm-sm0) + beta_vpd_gymno*(vpd-vpd0)
gymno_highsm <- data.frame(logrw, rw = exp(logrw), vpd, sm = high_sm, clade = "Gymnosperms")

predictions_vpd <- rbind(angio_lowsm, gymno_lowsm, angio_highsm, gymno_highsm)

ggplot() +
  annotate("rect", xmin = low_vpd-0.9, xmax = low_vpd+0.9, ymin = 0.83, ymax = 1.2, color = "#e8e832", fill = NA, linewidth = 1) +
  annotate("rect", xmin = high_vpd-0.9, xmax = high_vpd+0.9, ymin = 0.19, ymax = 0.34, color = "#b73131", fill = NA, linewidth = 1) +
  geom_line(aes(x = vpd, y = rw, group = paste0(sm, clade), color = clade, linetype = as.character(sm)), data = predictions_vpd) +
  geom_point(aes(x = vpd, y = rw, color = clade), data = predictions_vpd[predictions_vpd$vpd %in% c(low_vpd, high_vpd),]) +
  geom_boxplot(aes(x = vpd, y = 0), outliers = FALSE, width = 0.1, data = vpd_sm) +
  scale_color_manual(values = c('#27278f', '#278f27'), breaks = c("Gymnosperms", "Angiosperms"), name='Clade', guide = guide_legend(order = 2)) +
  scale_linetype_manual(values = c('solid', 'dashed'),labels = c('Low (19%)', 'High (32%)'), name = 'Soil\nmoisture', guide = guide_legend(order = 1)) +
  theme_classic() +
  labs(x = 'VPD (hPa)', y = 'Ring width (mm)')
  

low_vpd <- as.numeric(quantile(data$vpd_obs, 0.05)) # 3.67%
sm <- c(seq(range(data$sm_obs)[1], range(data$sm_obs)[2], 0.1), low_sm, high_sm)
logrw <- alpha + beta_sm_angio*(sm-sm0) + beta_vpd_angio*(low_vpd-vpd0)
angio_lowvpd <- data.frame(logrw, rw = exp(logrw), vpd = round(low_vpd,1), sm, clade = "Angiosperms")
logrw <- alpha + beta_sm_gymno*(sm-sm0) + beta_vpd_gymno*(low_vpd-vpd0)
gymno_lowvpd <- data.frame(logrw, rw = exp(logrw), vpd = round(low_vpd,1), sm, clade = "Gymnosperms")

high_vpd <- as.numeric(quantile(data$vpd_obs, 0.95)) # 18.01%
sm <- c(seq(range(data$sm_obs)[1], range(data$sm_obs)[2], 0.1), low_sm, high_sm)
logrw <- alpha + beta_sm_angio*(sm-sm0) + beta_vpd_angio*(high_vpd-vpd0)
angio_highvpd <- data.frame(logrw, rw = exp(logrw), vpd = round(high_vpd,1), sm, clade = "Angiosperms")
logrw <- alpha + beta_sm_gymno*(sm-sm0) + beta_vpd_gymno*(high_vpd-vpd0)
gymno_highvpd <- data.frame(logrw, rw = exp(logrw), vpd = round(high_vpd,1), sm , clade = "Gymnosperms")

predictions_sm <- rbind(angio_lowvpd, gymno_lowvpd, angio_highvpd, gymno_highvpd)

ggplot() +
  geom_line(aes(x = sm, y = rw, group = paste0(vpd, clade), color = clade, linetype = as.character(vpd)), data = predictions_sm) +
  geom_point(aes(x = sm, y = rw, color = clade), data = predictions_sm[predictions_sm$sm %in% c(low_sm, high_sm),]) +
  geom_boxplot(aes(x = sm, y = 0), outliers = FALSE, width = 0.1, data = vpd_sm) +
  scale_color_manual(values = c('#27278f', '#278f27'), breaks = c("Gymnosperms", "Angiosperms"), name='Clade', guide = guide_legend(order = 2)) +
  scale_linetype_manual(values = c('solid', 'dashed'), breaks = c(3.7, 18), labels = c('Low (3.6%)', 'High (18%)'), name = 'VPD', guide = guide_legend(order = 1)) +
  theme_classic() +
  labs(x = 'Soil moisture (m3.m-3)', y = 'Ring width (mm)') +
  annotate("rect", xmin = low_sm-0.9, xmax = low_sm+0.9, ymin = 0.13, ymax = 0.4, color = "#b73131", fill = NA, linewidth = 1) +
  annotate("rect", xmin = high_sm-0.9, xmax = high_sm+0.9, ymin = 0.83, ymax = 1.2, color = "#e8e832", fill = NA, linewidth = 1)



# What could happen with interactions?

beta_smvpd <- -0.01

low_sm <- as.numeric(quantile(data$sm_obs, 0.05)) # 19.10%
vpd <- c(seq(range(data$vpd_obs)[1], range(data$vpd_obs)[2], 0.1), low_vpd, high_vpd)
logrw <- alpha + beta_sm_angio*(low_sm-sm0) + beta_vpd_angio*(vpd-vpd0)+beta_smvpd*(vpd-vpd0)*(low_sm-sm0)
angio_lowsm <- data.frame(logrw, rw = exp(logrw), vpd, sm = low_sm, clade = "Angiosperms")
logrw <- alpha + beta_sm_gymno*(low_sm-sm0) + beta_vpd_gymno*(vpd-vpd0)+beta_smvpd*(vpd-vpd0)*(low_sm-sm0)
gymno_lowsm <- data.frame(logrw, rw = exp(logrw), vpd, sm = low_sm, clade = "Gymnosperms")

high_sm <- as.numeric(quantile(data$sm_obs, 0.95)) # 32.17%
vpd <- c(seq(range(data$vpd_obs)[1], range(data$vpd_obs)[2], 0.1), low_vpd, high_vpd)
logrw <- alpha + beta_sm_angio*(high_sm-sm0) + beta_vpd_angio*(vpd-vpd0)+beta_smvpd*(vpd-vpd0)*(high_sm-sm0)
angio_highsm <- data.frame(logrw, rw = exp(logrw), vpd, sm = high_sm, clade = "Angiosperms")
logrw <- alpha + beta_sm_gymno*(high_sm-sm0) + beta_vpd_gymno*(vpd-vpd0)+beta_smvpd*(vpd-vpd0)*(high_sm-sm0)
gymno_highsm <- data.frame(logrw, rw = exp(logrw), vpd, sm = high_sm, clade = "Gymnosperms")

predictions_vpd <- rbind(angio_lowsm, gymno_lowsm, angio_highsm, gymno_highsm)

ggplot() +
  annotate("rect", xmin = low_vpd-0.9, xmax = low_vpd+0.9, ymin = 1.2, ymax = 1.6, color = "#e8e832", fill = NA, linewidth = 1) +
  annotate("rect", xmin = high_vpd-0.9, xmax = high_vpd+0.9, ymin = 0.1, ymax = 0.4, color = "#b73131", fill = NA, linewidth = 1) +
  geom_line(aes(x = vpd, y = rw, group = paste0(sm, clade), color = clade, linetype = as.character(sm)), data = predictions_vpd) +
  geom_point(aes(x = vpd, y = rw, color = clade), data = predictions_vpd[predictions_vpd$vpd %in% c(low_vpd, high_vpd),]) +
  geom_boxplot(aes(x = vpd, y = 0), outliers = FALSE, width = 0.1, data = vpd_sm) +
  scale_color_manual(values = c('#27278f', '#278f27'), breaks = c("Gymnosperms", "Angiosperms"), name='Clade', guide = guide_legend(order = 2)) +
  scale_linetype_manual(values = c('solid', 'dashed'),labels = c('Low (19%)', 'High (32%)'), name = 'Soil\nmoisture', guide = guide_legend(order = 1)) +
  theme_classic() +
  labs(x = 'VPD (hPa)', y = 'Ring width (mm)')


low_vpd <- as.numeric(quantile(data$vpd_obs, 0.05)) # 3.67%
sm <- c(seq(range(data$sm_obs)[1], range(data$sm_obs)[2], 0.1), low_sm, high_sm)
logrw <- alpha + beta_sm_angio*(sm-sm0) + beta_vpd_angio*(low_vpd-vpd0)+beta_smvpd*(low_vpd-vpd0)*(sm-sm0)
angio_lowvpd <- data.frame(logrw, rw = exp(logrw), vpd = round(low_vpd,1), sm, clade = "Angiosperms")
logrw <- alpha + beta_sm_gymno*(sm-sm0) + beta_vpd_gymno*(low_vpd-vpd0)+beta_smvpd*(low_vpd-vpd0)*(sm-sm0)
gymno_lowvpd <- data.frame(logrw, rw = exp(logrw), vpd = round(low_vpd,1), sm, clade = "Gymnosperms")

high_vpd <- as.numeric(quantile(data$vpd_obs, 0.95)) # 18.01%
sm <- c(seq(range(data$sm_obs)[1], range(data$sm_obs)[2], 0.1), low_sm, high_sm)
logrw <- alpha + beta_sm_angio*(sm-sm0) + beta_vpd_angio*(high_vpd-vpd0)+beta_smvpd*(high_vpd-vpd0)*(sm-sm0)
angio_highvpd <- data.frame(logrw, rw = exp(logrw), vpd = round(high_vpd,1), sm, clade = "Angiosperms")
logrw <- alpha + beta_sm_gymno*(sm-sm0) + beta_vpd_gymno*(high_vpd-vpd0)+beta_smvpd*(high_vpd-vpd0)*(sm-sm0)
gymno_highvpd <- data.frame(logrw, rw = exp(logrw), vpd = round(high_vpd,1), sm , clade = "Gymnosperms")

predictions_sm <- rbind(angio_lowvpd, gymno_lowvpd, angio_highvpd, gymno_highvpd)

ggplot() +
  annotate("rect", xmin = low_sm-0.9, xmax = low_sm+0.9, ymin = 0.4, ymax = 0.6, color = "#b73131", fill = NA, linewidth = 1) +
  annotate("rect", xmin = high_sm-0.9, xmax = high_sm+0.9, ymin = 1.2, ymax = 1.6, color = "#e8e832", fill = NA, linewidth = 1) +
  geom_line(aes(x = sm, y = rw, group = paste0(vpd, clade), color = clade, linetype = as.character(vpd)), data = predictions_sm) +
  geom_point(aes(x = sm, y = rw, color = clade), data = predictions_sm[predictions_sm$sm %in% c(low_sm, high_sm),]) +
  geom_boxplot(aes(x = sm, y = 0), outliers = FALSE, width = 0.2, data = vpd_sm) +
  scale_color_manual(values = c('#27278f', '#278f27'), breaks = c("Gymnosperms", "Angiosperms"), name='Clade', guide = guide_legend(order = 2)) +
  scale_linetype_manual(values = c('solid', 'dashed'), breaks = c(3.7, 18), labels = c('Low (3.6%)', 'High (18%)'), name = 'VPD', guide = guide_legend(order = 1)) +
  theme_classic() +
  labs(x = 'Soil moisture (m3.m-3)', y = 'Ring width (mm)')
  



