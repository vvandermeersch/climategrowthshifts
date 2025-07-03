
library(ggplot2)
library(patchwork)
library(ggtree)


probs <- c(0.1, 0.5, 0.9)
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_gdd[', sp, ']'))
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles <- sapply(1:length(names), calc)
betagdd_est <- data.frame(uniq_species_ids, t(quantiles))
names(betagdd_est)[1:4] <- c('shortname', paste0('beta_gdd_', c('Q10', 'Q50', 'Q90')))

probs <- c(0.1, 0.5, 0.9)
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_sm[', sp, ']'))
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles <- sapply(1:length(names), calc)
betasm_est <- data.frame(uniq_species_ids, t(quantiles))
names(betasm_est)[1:4] <- c('shortname', paste0('beta_sm_', c('Q10', 'Q50', 'Q90')))

probs <- c(0.1, 0.5, 0.9)
names <- sapply(1:data$N_species,
                function(sp) paste0('beta_vpd[', sp, ']'))
calc <- function(n) {
  util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
}
quantiles <- sapply(1:length(names), calc)
betavpd_est <- data.frame(uniq_species_ids, t(quantiles))
names(betavpd_est)[1:4] <- c('shortname', paste0('beta_vpd_', c('Q10', 'Q50', 'Q90')))

beta_estimates <- merge(merge(betagdd_est, betasm_est), betavpd_est)

beta_estimates <- merge(beta_estimates, sppfull[,c('shortname', 'phylo.name')])
names(beta_estimates)[names(beta_estimates)=='phylo.name'] <- 'label'
tree <- ggtree::ggtree(phy.plants.here) 
datatree <- tree$data
datatree <- merge(datatree, beta_estimates, all.x = TRUE)

kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF",
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")

tree <- ggtree(datatree) +
  geom_tiplab(color = 'grey30') +
  ggplot2::scale_color_gradientn(colors = kippenberger) +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, nrow(beta_estimates)+1), expand = FALSE)

gdd <- ggplot(data = datatree) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey80') +
  geom_pointrange(aes(x = beta_gdd_Q50, xmin = beta_gdd_Q10, xmax= beta_gdd_Q90,
                      color = beta_gdd_Q50, y = y),
                  size = 0.2) +
  scale_color_gradientn(colors = kippenberger, limits = c(-0.2,0.2)) +
  guides(
    color = guide_colorbar(order = 0,
                           frame.colour = "grey20", ticks.colour = NA,
                           frame.linewidth = 0.2)) +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
        legend.position = "inside", legend.text = element_text(size = 8),
        legend.position.inside = c(0.85,0.15), legend.title = element_blank(),
        legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
  coord_cartesian(xlim = c(-0.5,0.5), ylim = c(0, nrow(beta_estimates)+1), expand = FALSE) +
  labs(x = 'GDD')

sm <- ggplot(data = datatree) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey80') +
  geom_pointrange(aes(x = beta_sm_Q50, xmin = beta_sm_Q10, xmax= beta_sm_Q90,
                      color = beta_sm_Q50, y = y),
                  size = 0.2) +
  scale_color_gradientn(colors = kippenberger, limits = c(-0.5,0.5)) +
  guides(
    color = guide_colorbar(order = 0,
                           frame.colour = "grey20", ticks.colour = NA,
                           frame.linewidth = 0.2)) +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
        legend.position = "inside", legend.text = element_text(size = 8),
        legend.position.inside = c(0.85,0.15), legend.title = element_blank(),
        legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
  coord_cartesian(xlim = c(-0.5,0.5), ylim = c(0, nrow(beta_estimates)+1), expand = FALSE) +
  labs(x = 'Soil moisture')

vpd <- ggplot(data = datatree) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey80') +
  geom_pointrange(aes(x = beta_vpd_Q50, xmin = beta_vpd_Q10, xmax= beta_vpd_Q90, 
                      color = beta_vpd_Q50, y = y),
                  size = 0.2) +
  scale_color_gradientn(colors = kippenberger, limits = c(-0.2,0.2)) +
  guides(
    color = guide_colorbar(order = 0,
                           frame.colour = "grey20", ticks.colour = NA,
                           frame.linewidth = 0.2)) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = margin(t = 5.5, b = 5.5, r = 12, l = 10),
        legend.position = "inside", legend.text = element_text(size = 8),
        legend.position.inside = c(0.85,0.15), legend.title = element_blank(),
        legend.direction="vertical", legend.key.width  = unit(8, "pt")) +
  coord_cartesian(xlim = c(-0.5,0.5), ylim = c(0, nrow(beta_estimates)+1), expand = FALSE) +
  labs(x = 'VPD')
  

slopes_estimates <- tree + gdd + sm + vpd + plot_layout(widths = c(1, 0.7, 0.7, 0.7), nrow = 1)

ggsave(slopes_estimates, file = file.path(wd, 'figures', 'model28june_nopooling', 'slopes_estimates.pdf'), 
       width = 297, height = 210, units = "mm")


