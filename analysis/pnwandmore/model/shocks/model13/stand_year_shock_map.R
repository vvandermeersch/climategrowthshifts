for(s in 1:data$N_stands){
  print(s)
  
  for(y in stand_years[[s]]){
    
    propname <- paste0('prop_sck_conc_states[', s, ',', y, ']')
    gq_samples[[propname]] <- 0
    
    obs <- unlist(sapply(which(data$stand_idxs == s), function(i) data$tree_start_idxs[i]:data$tree_end_idxs[i]))
    obs_y <- which(data$all_years_idxs == y)
    obs_y <- intersect(obs, obs_y)
    
    for(i in obs_y){
      name_i <- paste0('sck_conc_state[', i, ']')
      gq_samples[[propname]] <- gq_samples[[propname]] + gq_samples[[name_i]]
    }
    gq_samples[[propname]] <- gq_samples[[propname]]/length(obs_y)
    
  }
  
}

library(terra)
northamerica <- geodata::gadm(country = geodata::country_codes("North America")$ISO3, level = 1, resolution = 2, path = file.path(wd, 'data'))


par(mfrow = c(2,3), mar = c(2.5,3,2,1), cex.main = 1, cex.lab = 1.2, cex.main = 2)
for(y in c(1956,1977,1992,2002,2004,2018)){
  
  data_plot <- data.frame(x = data$uniq_stand_lon, y = data$uniq_stand_lat,
                          prop_sck = sapply(1:data$N_stands, function(s){
                            if(is.null(gq_samples[[paste0('prop_sck_conc_states[', s, ',', which(data$all_years == y), ']')]])){return(NA)}
                            util$ensemble_mcmc_quantile_est(gq_samples[[paste0('prop_sck_conc_states[', s, ',', which(data$all_years == y), ']')]], 0.5)}))
  
  
  pal <- colorRampPalette(c("#FFEC9DFF", "#FAC881FF", "#F4A464FF", "#E87444FF", "#D9402AFF", "#BF2729FF", "#912534FF", "#64243EFF"))
  cols <- pal(100)
  z <- data_plot$prop_sck
  z_col <- cols[cut(z, breaks = seq(0,1,length =100), include.lowest = TRUE)]
  z_col<- ifelse(z == 0, "white", z_col)
  z_col<- ifelse(is.na(z), "grey90", z_col)
  
  
  
  plot(
    # asp = 1,
    x = NULL, y = NULL,
    type = "n",
    xlim = c(-127, -102),
    ylim = c(30, 50),
    xaxs = "i",  yaxs = "i",
    xlab = "", ylab= "",
    main = y
  )
  lines(crop(northamerica, ext(c(-127, -102, 30,50))), col = 'grey50')
  points(x = data$uniq_stand_lon, y = data$uniq_stand_lat, pch = 20, col = "grey20", cex = 2)
  points(x = data$uniq_stand_lon, y = data$uniq_stand_lat, pch = 20, col = z_col, cex = 1.6)
  points(x = data$uniq_stand_lon, y = data$uniq_stand_lat, pch = 4, col = ifelse(is.na(z), "grey20", NA), cex = 1)
  
  
  
  # legend position
  x0 <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  y0 <- par("usr")[3]
  legend_colors <- pal(10)
  legend_values <- round(seq(0, 1, length.out = 11), 2)
  rect(
    xleft = -124.5, xright = -124,
    ybottom = 31 + seq(0, 4, length.out = 10),
    ytop    = 31 + seq(1, 5, length.out = 10),
    col = legend_colors,
    border = NA
  )
  text(
    x = -124+0.2,
    y = 31 + seq(0, 5, length.out = 11),
    labels = legend_values,
    adj = 0, cex = 0.7
  )
  
}



