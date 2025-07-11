rm(list = ls())
wd <- "~/projects/climategrowthshifts/analysis/pnw"
library(ggplot2)
setwd(file.path(wd, 'model'))
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
setwd(wd)

# Model data
data <- readRDS(file.path(wd, 'output/model', 'data_28june2025.rds'))

# Posterior quantification - diagnostics
fit <- readRDS(file.path(wd, 'output/model', 'fit_28june2025_nopooling.rds')) 


util$plot_realizations(samples, names, data$all_years,N_plots=50,
                       xlab="Year",
                       display_ylim=c(-3, 3), display_xlim = c(1980, 2015))
samples <- util$extract_expectand_vals(fit)


stands <- 1:data$N_stands
fshs_df <- data.frame()

for(stand in stands){
  
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  names <- util$check_expectand_names(names, samples)
  
  # Extract function values
  
  fs <- t(sapply(names, function(name)
    c(samples[[name]], recursive=TRUE)))
  N <- dim(fs)[1]
  I <- dim(fs)[2]
  J <- 20
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    fshs_df <- rbind(fshs_df, data.frame(stand = stand, real = j, year = data$all_years, fsh = f))
  }
}




ggplot(data = fshs_df[fshs_df$stand == 129,]) +
  geom_line(aes(x = year, y = fsh, color = stand, group = paste0(stand,'.',real)),
            linewidth = 0.3, alpha = 0.8) +
  annotate(geom = 'rect', xmin = 2013, xmax = max(data$all_years), ymin = -6.2, ymax = 4, fill = 'white', alpha = 0.7) +
  geom_vline(aes(xintercept = 2002), color = "white", linewidth = 1) +
  geom_vline(aes(xintercept = 2002), linetype = 'dashed') +
  # stat_summary(aes(x = year, y = fsh, group=1), fun=median, colour="white", 
  #              geom="line",group=1, linewidth = 2) +
  # stat_summary(aes(x = year, y = fsh, group=1), fun=median, colour="darkred", 
  #              geom="line",group=1, linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  theme(legend.position = 'none')



ggplot(data = fshs_df) +
  geom_line(aes(x = year, y = fsh, color = stand, group = paste0(stand,'.',real)),
            linewidth = 0.3, alpha = 0.8) +
  annotate(geom = 'rect', xmin = 2013, xmax = max(data$all_years), ymin = -6.2, ymax = 4, fill = 'white', alpha = 0.7) +
  geom_vline(aes(xintercept = 1986), color = "white", linewidth = 1) +
  geom_vline(aes(xintercept = 1986), linetype = 'dashed') +
  geom_vline(aes(xintercept = 1991), color = "white", linewidth = 1) +
  geom_vline(aes(xintercept = 1991), linetype = 'dashed') +
  stat_summary(aes(x = year, y = fsh, group=1), fun=median, colour="white",
               geom="line",group=1, linewidth = 2) +
  stat_summary(aes(x = year, y = fsh, group=1), fun=median, colour="darkred",
               geom="line",group=1, linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  theme(legend.position = 'none')

ggplot(data = pred) +
  geom_line(aes(x = year, y = vpd_mjj, color = ID, group = ID),
            linewidth = 0.3, alpha = 0.3) +
  # geom_vline(aes(xintercept = 1986), color = "white", linewidth = 1) +
  # geom_vline(aes(xintercept = 1986), linetype = 'dashed') +
  # geom_vline(aes(xintercept = 1991), color = "white", linewidth = 1) +
  # geom_vline(aes(xintercept = 1991), linetype = 'dashed') +
  geom_vline(aes(xintercept = 2002), color = "white", linewidth = 1) +
  geom_vline(aes(xintercept = 2002), linetype = 'dashed') +
  stat_summary(aes(x = year, y = vpd_mjj, group=1), fun=median, colour="white",
               geom="line",group=1, linewidth = 2) +
  stat_summary(aes(x = year, y = vpd_mjj, group=1), fun=median, colour="darkred",
               geom="line",group=1, linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  theme(legend.position = 'none')

ggplot(data = pred) +
  geom_line(aes(x = year, y = soilmoist_mjj*100, color = ID, group = ID),
            linewidth = 0.3, alpha = 0.3) +
  # geom_vline(aes(xintercept = 1986), color = "white", linewidth = 1) +
  # geom_vline(aes(xintercept = 1986), linetype = 'dashed') +
  # geom_vline(aes(xintercept = 1991), color = "white", linewidth = 1) +
  # geom_vline(aes(xintercept = 1991), linetype = 'dashed') +
  geom_vline(aes(xintercept = 2002), color = "white", linewidth = 1) +
  geom_vline(aes(xintercept = 2002), linetype = 'dashed') +
  stat_summary(aes(x = year, y = soilmoist_mjj*100, group=1), fun=median, colour="white",
               geom="line",group=1, linewidth = 2) +
  stat_summary(aes(x = year, y = soilmoist_mjj*100, group=1), fun=median, colour="darkred",
               geom="line",group=1, linewidth = 1) +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  theme(legend.position = 'none')
            