

util$plot_realizations(samples, names, data$all_years,N_plots=50,
                       xlab="Year",
                       display_ylim=c(-3, 3), display_xlim = c(1980, 2015))

par(mar = c(2,2,1,1))
plot(1, type="n", xlim = c(1980, 2000), ylim = c(-5,5))
for(stand in 1:10){
  standname <- uniq_stand_ids[stand]
  lat <- unique(datasets[datasets$grouped_stand == standname, "north_lat"])
  
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  names <- util$check_expectand_names(names, samples)
  
  # Extract function values
  fs <- t(sapply(names, function(name)
    c(samples[[name]], recursive=TRUE)))
  N <- dim(fs)[1]
  I <- dim(fs)[2]
  J <- 1
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  
  plot_xs <- data$all_years
  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    lines(plot_xs, f, lwd=1)
  }
}


latstand <- datasets[,c("grouped_stand","north_lat", "east_lon")]
latstand$north_lat <- round(latstand$north_lat,2)
latstand$east_lon <- round(latstand$east_lon,2)
latstand <- unique(latstand)
latstand <- latstand[latstand$grouped_stand %in% unique(raw_data$grouped_stand),]
latstand$stand_number <- 1:nrow(latstand)

stands48lat <- latstand[latstand$north_lat >= 48,]
stands35lat <- latstand[latstand$north_lat <= 35,]

par(mfrow = c(2,1),mar = c(2,2,1,1))
plot(1, type="n", xlim = c(1980, 2000), ylim = c(-3,3), , main = "Latitude >= 48deg")
for(stand in stands48lat$stand_number){
  
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  names <- util$check_expectand_names(names, samples)
  
  # Extract function values
  fs <- t(sapply(names, function(name)
    c(samples[[name]], recursive=TRUE)))
  N <- dim(fs)[1]
  I <- dim(fs)[2]
  J <- 1
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  
  plot_xs <- data$all_years
  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    lines(plot_xs, f, lwd=1, co = "darkblue")
  }
}
plot(1, type="n", xlim = c(1980, 2000), ylim = c(-3,3), main = "Latitude <= 33deg")

for(stand in stands33lat$stand_number){
  
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  names <- util$check_expectand_names(names, samples)
  
  # Extract function values
  fs <- t(sapply(names, function(name)
    c(samples[[name]], recursive=TRUE)))
  N <- dim(fs)[1]
  I <- dim(fs)[2]
  J <- 1
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  
  plot_xs <- data$all_years
  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    lines(plot_xs, f, lwd=1, co = "darkred")
  }
}


library(rnaturalearth)
library(terra)
library(tidyterra)

world <- vect(ne_countries(scale = "medium", returnclass = "sf"))
ggplot() +
  geom_spatvector(data = world, color = 'white', fill = 'grey90') +
  geom_point(data = latstand, aes(x = east_lon, y = north_lat, color = species_code), col = "black", size = 0.3) +
  geom_point(data = stands48lat, aes(x = east_lon, y = north_lat, color = species_code), col = "darkblue") +
  geom_point(data = stands33lat, aes(x = east_lon, y = north_lat, color = species_code), col = "darkred") +
  ggplot2::theme_minimal() +
  theme(axis.title = ggplot2::element_blank()) +
  ggplot2::coord_sf(xlim = c(-127, -100), ylim = c(30,52))



stands_sierranevada <- latstand[latstand$north_lat <= 40 & latstand$north_lat> 35  & latstand$east_lon <= -116 ,]

ggplot() +
  geom_spatvector(data = world, color = 'white', fill = 'grey90') +
  geom_point(data = stands_sierranevada, aes(x = east_lon, y = north_lat, color = species_code), col = "black", size = 0.3) +
  ggplot2::coord_sf(xlim = c(-127, -100), ylim = c(30,52))


par(mfrow = c(1,1),mar = c(2,2,1,1))
plot(1, type="n", xlim = c(1980, 2002), ylim = c(-3,3), main = "Latitude <= 33deg")
fshorts <- data.frame()
for(stand in stands_sierranevada$stand_number){
  
  names <- sapply(1:length(data$all_years),
                  function(n) paste0('f_sh[', stand,',', n, ']'))
  names <- util$check_expectand_names(names, samples)
  
  # Extract function values
  fs <- t(sapply(names, function(name)
    c(samples[[name]], recursive=TRUE)))
  N <- dim(fs)[1]
  I <- dim(fs)[2]
  J <- 1
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  
  plot_xs <- data$all_years
  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    lines(plot_xs, f, lwd=1, co = "darkred")
  }
  
  fshorts <- rbind(fshorts, data.frame(f, stand, year = plot_xs))
}

clim_sierranevada <- clim_pred[clim_pred$dataset %in%datasets[datasets$grouped_stand %in% stands_sierranevada$grouped_stand, 'dataset'],]


ggplot() +
  # ggplot2::geom_line(data = clim_sierranevada, aes(x = year, y = vpd_mjj, group = dataset), alpha = 0.1) +
  ggplot2::stat_summary(data = clim_sierranevada, aes(y = vpd_mjj, x = year, group=1), fun=mean, colour="darkred", geom="line",group=1) +
  ggplot2::stat_summary(data = clim_sierranevada, aes(y = soilmoist_mjj, x = year, group=1), fun=mean, colour="darkblue", geom="line",group=1) +
  ggplot2::geom_line(data = fshorts, aes(x = year, y = f*5+15, group = stand), color = "black", alpha = 0.3) +
  ggplot2::theme_classic() +
  ggplot2::xlim(c(1980,2002)) +
  ggplot2::geom_vline(xintercept = 1995, color = 'darkgrey') 

ggplot(data = clim_sierranevada) +
  ggplot2::geom_line(aes(x = year, y = soilmoist_mjj, group = dataset))

ggplot(data = clim_sierranevada) +
  ggplot2::geom_line(aes(x = year, y = gdd, group = dataset))
  