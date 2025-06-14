rm(list = ls())

library(terra)
library(ggplot2)

kippenberger <- c("#8B174DFF", "#AE2565FF", "#C1447EFF", "#D06C9BFF", "#DA9FB8FF", 
                  "#ADBE7CFF", "#8BA749FF", "#6E8537FF", "#4F5F28FF", "#343D1FFF")

wd <- "/home/victor/projects/climategrowthshifts/analysis/mountrainier/data"

plots <- read.csv(file.path(wd, "treerings_ailene", "tree_plot_climate_temp.csv"))
plots <- vect(unique(plots[,c("Longitude", "Latitude", "Plot", "Elevation")]), geom=c("Longitude", "Latitude"))

plotsID <- data.frame(ID = 1:16, plotname = plots$Plot, alt = plots$Elevation)

years <- seq(1980,2023,1)
vars <- list("tair" = "Tair_f_tavg", 
             "swe" = "SWE_tavg",
             "soilmoist_010" = "SoilMoi00_10cm_tavg",
             "soilmoist_1040" = "SoilMoi10_40cm_tavg",
             "soilmoist_40100" = "SoilMoi40_100cm_tavg",
             "soilmoist_100200" = "SoilMoi100_200cm_tavg")
datdf <- data.frame()

for(year in years){
  
  cat(paste0(year,"\n"))
  # files <- list.files(path = , pattern = as.character(year), full.names = TRUE)
  
  dates <- format(seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")),1), "%Y%m%d")
  files <- file.path(wd, "climate","wldas", paste0("WLDAS_NOAHMP001_DA1_",dates,".D10.nc.SUB.nc4"))
  
  
  # if(length(files) < 365){stop("Problem with no. of files")}
  if(any(!file.exists(files))){stop("Problem with no. of files")}
  dat <- rast(files)
  
  daty <- data.frame()
  
  for(v in vars){
    newv <- names(vars)[which(vars==v)]
    cat(paste0("   ", newv,"\n"))
    
    datv <- extract(dat[v], plots)
    datv <- as.data.frame(tidyr::pivot_longer(datv, cols = -ID, names_to = "var", values_to = "value"))
    
    if(newv == "tair"){
      datv$value <- datv$value  -273.15 # K to C
    }
    
    datv$var <- newv
    daty <- rbind(daty, datv)
  }
  
  daty$date <- rep(seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")),1), length(unique(daty$ID))*length(vars))
  daty$doy <- lubridate::yday(daty$date)
  datdf <- rbind(datdf, daty)
  
}

# reshape(datdf, direction = "wide", timevar = "date", idvar = "ID")

ggplot(data = datdf[datdf$var == "swe" & lubridate::year(datdf$date) %in% c(1980:2023) ,]) +
  geom_line(aes(x = date, y = value, group = ID, color = ID),
            linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) + 
  theme_bw() +
  scale_x_date(expand = c(0,0)) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = "DOY", y = "SWE (mm)")

swe <- extract(yearnc["SWE_tavg"], plots)
swe <- as.data.frame(tidyr::pivot_longer(swe, cols = -ID, names_to = "var", values_to = "swe"))
swe$doy <- rep(3:365, length(unique(swe$ID)))

ggplot(data = swe) +
  geom_line(aes(x = doy, y = swe, group = ID, color = ID),
            linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) + 
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = "DOY", y = "SWE (mm)")

soilmoist010 <- extract(yearnc["SoilMoi00_10cm_tavg"], plots)
soilmoist010 <- as.data.frame(tidyr::pivot_longer(soilmoist010, cols = -ID, names_to = "var", values_to = "soilmoist010"))
soilmoist010$doy <- rep(3:365, length(unique(soilmoist010$ID)))

ggplot(data = soilmoist010) +
  geom_line(aes(x = doy, y = soilmoist010, group = ID, color = ID),
            linewidth = 0.3) +
  scale_color_gradientn(colors = kippenberger) + 
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = 'none',
        panel.grid = element_blank()) +
  labs(x = "DOY", y = "Soil moisture, 0-10cm (m3.m-3)")




