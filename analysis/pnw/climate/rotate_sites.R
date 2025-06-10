

for(yr in 1980:2024){
  sites <- vect(df_pheno[df_pheno$year == yr,], geom = c("longitude", "latitude"), crs = "EPSG:4326")
  sites <- shift(rotate(sites, longitude = 0),360) # cannot use this, as I don't have the last version of terra
  
  saveRDS(sites, file.path(wd, 'output', 'pheno', paste0("sitesrotated_",yr,".rds")))
}

