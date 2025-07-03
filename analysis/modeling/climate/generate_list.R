

dates <- seq(as.Date("1979-01-03"), as.Date("2023-12-31"), by = "day")
dates <- seq(as.Date("1980-01-01"), as.Date("2023-12-31"), by = "day")

dates <- c(as.Date("1992-12-07"), as.Date("1992-12-08"), as.Date("1992-12-09"))

dates <- seq(as.Date("1980-01-01"), as.Date("2023-12-31"), by = "day")

listdates <- data.frame()
for(d in dates){
  d <- as.Date(d)
  file <- paste0("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FWLDAS%2FWLDAS_NOAHMP001_DA1.D1.0%2F",format(d, '%Y'),"%2F",format(d, '%m'),"%2FWLDAS_NOAHMP001_DA1_",format(d, '%Y%m%d'),".D10.nc&BBOX=25.065%2C-124.925%2C52.925%2C-89.025&SERVICE=L34RS_LDAS&SHORTNAME=WLDAS_NOAHMP001_DA1&LABEL=WLDAS_NOAHMP001_DA1_",format(d, '%Y%m%d'),".D10.nc.SUB.nc4&VARIABLES=SoilMoi00_10cm_tavg%2CSoilMoi10_40cm_tavg%2CSoilMoi100_200cm_tavg%2CSoilMoi40_100cm_tavg%2CTair_f_tavg&DATASET_VERSION=D1.0&FORMAT=bmM0Lw&VERSION=1.02")
  if(format(d, '%m%d') == "0101"){
    file <- paste0("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FWLDAS%2FWLDAS_NOAHMP001_DA1.D1.0%2F",as.numeric(format(d, '%Y'))-1,"%2F",12,"%2FWLDAS_NOAHMP001_DA1_",format(d, '%Y%m%d'),".D10.nc&BBOX=25.065%2C-124.925%2C52.925%2C-89.025&SERVICE=L34RS_LDAS&SHORTNAME=WLDAS_NOAHMP001_DA1&LABEL=WLDAS_NOAHMP001_DA1_",format(d, '%Y%m%d'),".D10.nc.SUB.nc4&VARIABLES=SoilMoi00_10cm_tavg%2CSoilMoi10_40cm_tavg%2CSoilMoi100_200cm_tavg%2CSoilMoi40_100cm_tavg%2CTair_f_tavg&DATASET_VERSION=D1.0&FORMAT=bmM0Lw&VERSION=1.02")
  }else if(format(d, '%d') == "01"){
    file <- paste0("https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FWLDAS%2FWLDAS_NOAHMP001_DA1.D1.0%2F",format(d, '%Y'),"%2F",format(d %m-% months(1), '%m'),"%2FWLDAS_NOAHMP001_DA1_",format(d, '%Y%m%d'),".D10.nc&BBOX=25.065%2C-124.925%2C52.925%2C-89.025&SERVICE=L34RS_LDAS&SHORTNAME=WLDAS_NOAHMP001_DA1&LABEL=WLDAS_NOAHMP001_DA1_",format(d, '%Y%m%d'),".D10.nc.SUB.nc4&VARIABLES=SoilMoi00_10cm_tavg%2CSoilMoi10_40cm_tavg%2CSoilMoi100_200cm_tavg%2CSoilMoi40_100cm_tavg%2CTair_f_tavg&DATASET_VERSION=D1.0&FORMAT=bmM0Lw&VERSION=1.02")
  }
  
  listdates <- rbind(listdates, file)
}
write.table(listdates, file = "/home/victor/projects/climategrowthshifts/analysis/modeling/data/climate/list2.txt", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --content-disposition -i 'list3.txt'


