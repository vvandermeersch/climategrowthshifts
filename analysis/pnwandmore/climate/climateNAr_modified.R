# I have to modify the function so it can work
ClimateNAr <- function(inputFile='Test/test.csv',periodList='Normal_1961_1990.nrm',varList=c('MAT','MAP'),outDir='c:/temp/',
                       prefix = 'clim'){
  ## Initial steps ================================================>
  # library(terra); library(data.table)
  
  #version verification =================--------------
  tryCatch({
    url <- "https://climatena.ca/downloads/ClimateNAr_version310.csv" #by changing the filename, the verification will fail.
    dat0 <- read.csv(url)
    message(dat0[1])
  }, warning = function(w) {
    stop("Version verification failed.")
  },  error = function(e) {
    stop("Version verification failed.")
  }) # end of tryCatch
  
  
  ## input coord file ***********************************=====================-------------------------------
  # determine inputFile type: CSV vs. raster
  # message('inputFile')
  #if(grepl('.tif|.asc|.TIF|.ASC', inputFile)){
  # inputFile_class='Raster'
  #}
  # if(grepl('.csv|.CSV', inputFile)){
  #inputFile_class='CSV'
  #}
  #inputFile_class
  #if(is.data.frame(inputFile)){inputFile_class='dataFrame'};inputFile_class
  
  if(is.character(inputFile)){
    if(grepl('.tif|.asc|.TIF|.ASC',inputFile)){inputFile_class='Raster'}else{inputFile_class='CSV'};inputFile_class
    if(grepl('.csv|.CSV',inputFile)){inputFile_class='CSV'};inputFile_class
  }
  if(is.data.frame(inputFile)){inputFile_class='dataFrame'};inputFile_class
  
  
  # CSV files
  if(inputFile_class=='CSV'){
    input <- as.data.table(readr::read_csv(inputFile,show_col_types = FALSE));head(input);dim(input)
    elev0 <- input[[5]];head(elev0)
    #xy_input0 <- input[,c(4,3)];names(xy_input0) <- c('x','y');head(xy_input0)
    xyl_input <- cbind(input[,c(4,3)],elev0);head(xyl_input)
    if(dim(input)[1]>1 & dim(input)[2]==5){input_valid=T}; input_valid
    if(dim(input)[1]>1000000){message('Your input file might be too big')}
  }
  
  # dataFrame
  if(inputFile_class=='dataFrame'){
    input <- as.data.table(inputFile)
    elev0 <- input[[5]];head(elev0)
    #xy_input0 <- input[,c(4,3)];names(xy_input0) <- c('x','y');head(xy_input0)
    xyl_input <- cbind(input[,c(4,3)],elev0)
    # periodList = periodList[1]
    if(dim(input)[1]>1 & dim(input)[2]==5){input_valid=T}; input_valid
    if(dim(input)[1]>1000000){message('Your input file might be too big')}
  }
  
  # raster file
  if(inputFile_class=='Raster'){
    input_tif <- terra::rast(inputFile);#input_tif #plot(input_tif)
    ext_input <- ext(input_tif);ext_input #for final raster output
    origin_input <- origin(input_tif);origin_input #for final raster output
    names(input_tif) <- 'elev0'
    xyl_input <- as.data.table(cbind(crds(input_tif,na.rm=T), values(input_tif,na.rm=T)));head(xyl_input);dim(xyl_input)
    if(dim(input_tif)[1]>10){input_valid=T}; input_valid
    if(dim(xyl_input)[1]>5000000){message('Your input file might be too big')}
    rm(input_tif);gc()
  }
  
  
  ## functions ====================
  # load system data ------
  load_object <- function(file) {
    tmp <- new.env()
    load(file = file, envir = tmp)
    tmp[[ls(tmp)[1]]]
  }
  
  # RH ---------
  calculate_RH <- function(Tmax, Tmin) {
    if (ncol(Tmax) != ncol(Tmin)) {
      stop("Tmax and Tmin must have the same number of columns")
    }
    n_months <- ncol(Tmax)
    n_rows <- nrow(Tmax)
    names_col <- gsub('Tmax','RH', names(Tmax))
    SVP_Tmax <- 0.6105 * exp((17.273 * Tmax) / (Tmax + 237.3))
    SVP_Tmin <- 0.6105 * exp((17.273 * Tmin) / (Tmin + 237.3))
    RH <- matrix(nrow = n_rows, ncol = n_months) # Initialize results matrix
    for(i in 1:n_months) {
      es_Tmax <- ifelse(Tmax[[i]] >= 0, SVP_Tmax[[i]], SVP_Tmax[[i]] * (1 + (Tmax[[i]] * 0.01)))
      es_Tmin <- ifelse(Tmin[[i]] >= 0, SVP_Tmin[[i]], SVP_Tmin[[i]] * (1 + (Tmin[[i]] * 0.01)))
      SVP_Tave <- (es_Tmax + es_Tmin) / 2
      RH[,i] <- 100 * es_Tmin / SVP_Tave
    }
    #colnames(RH) <- paste0('RH', sprintf("%02d", 1:n_months))
    colnames(RH) <- names_col
    RH <- data.table::as.data.table(round(RH,0))
    return(RH)
  }
  #rh4 <- calculate_RH(Tmax3, Tmin3);head(rh4)
  #rh4 <- calculate_RH(Tmax12, Tmin12);head(rh4)
  
  ## Eref, ErefR, CMD -------------------------------------------------------------
  #Tmax3=Tmax2; Tmin3=Tmin2; Tave3=Tave2;ppt3=ppt2
  calculate_Eref_CMD <- function(xyl_set, Tmax3, Tmin3, Tave3, ppt3) {
    # Input validation
    if(nrow(Tmax3) != nrow(Tmin3)) {
      stop("Tmax and Tmin must have the same number of rows")
    }
    n_rows <- nrow(Tmax3)
    n_cols <- ncol(Tmax3)
    names_Eref <- gsub('Tmax','Eref',names(Tmax3))
    names_ErefR <- gsub('Tmax','ErefR',names(Tmax3));names_ErefR
    names_CMD <- gsub('Tmax','CMD',names(Tmax3))
    digt_month <- as.numeric(gsub("\\D", "", names(Tmax3)))
    # Constants
    dsim <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    srdis <- c(1.032, 1.023, 1.01, 0.992, 0.977, 0.968, 0.968, 0.976, 0.991, 1.008, 1.023, 1.032)
    decln <- c(-0.3681, -0.2259, -0.0422, 0.1706, 0.332, 0.4076, 0.373, 0.2355, 0.0327, -0.1729, -0.3374, -0.4071)
    # Initialize output matrices
    Eref <- matrix(0, nrow = n_rows, ncol = n_cols)
    CMD <- matrix(0, nrow = n_rows, ncol = n_cols)
    ErefR <- matrix(0, nrow = n_rows, ncol = n_cols)
    
    # Latitude calculations
    lat_rad = xyl_set$y * 0.01745
    SinLat = sin(lat_rad)
    Coslat = cos(lat_rad)
    TanLat = tan(lat_rad)
    Lat_correction = (1.19 - 0.0067 * xyl_set$y)
    # Main calculations for each month
    if(n_cols==12){
      for(i in 1:12) {
        Scratch = (TanLat * tan(decln[i]))
        suppressWarnings({
          SunsetHrAngle = ifelse(Scratch >= -1 & Scratch <= 1, atan(Scratch / sqrt(Scratch * (-Scratch) + 1)) + 2 * atan(1), 3.14)
        })
        SunsetHrAngle = ifelse(Scratch < -1, 0, SunsetHrAngle)
        SolarTop_mm = 15.392 * srdis[i] * (SunsetHrAngle[[1]] * SinLat * sin(decln[i]) + (Coslat * cos(decln[i]) * sin(SunsetHrAngle[[1]])))
        Trange = Tmax3[[i]] - Tmin3[[i]]
        em = 0.0023 * dsim[i] * SolarTop_mm * (Tave3[[i]] + 17.8) * Trange ^ 0.5
        erefm0_values = em * Lat_correction
        erefr0_values = erefm0_values
        erefm0_values = ifelse(Tave3[[i]] < 0, 0, erefm0_values)
        cmd_values = ifelse(erefm0_values > ppt3[[i]], erefm0_values - ppt3[[i]], 0)
        # Store results in matrices
        Eref[,i] <- round(erefm0_values,0)
        CMD[,i] <- round(cmd_values,0)
        ErefR[,i] <- round(erefr0_values,0)
      }
    } else{
      i=12
      Scratch = (TanLat * tan(decln[i]))
      suppressWarnings({
        SunsetHrAngle = ifelse(Scratch >= -1 & Scratch <= 1, atan(Scratch / sqrt(Scratch * (-Scratch) + 1)) + 2 * atan(1), 3.14)
      })
      SunsetHrAngle = ifelse(Scratch < -1, 0, SunsetHrAngle)
      SolarTop_mm = 15.392 * srdis[i] * (SunsetHrAngle[[1]] * SinLat * sin(decln[i]) + (Coslat * cos(decln[i]) * sin(SunsetHrAngle[[1]])))
      Trange = Tmax3 - Tmin3
      em = 0.0023 * dsim[i] * SolarTop_mm * (Tave3 + 17.8) * Trange ^ 0.5
      erefm0_values = unlist(em * Lat_correction)
      erefr0_values <-  erefm0_values
      erefm0_values = unlist(ifelse(Tave3[[1]] < 0, 0, erefm0_values))
      cmd_values = ifelse(erefm0_values > ppt3[[1]], erefm0_values - ppt3[[1]], 0)
      # Store results in matrices
      Eref[,1] <- round(erefm0_values,0)
      CMD[,1] <- round(cmd_values,0)
      ErefR[,1] <- (round(erefr0_values,0))
    }
    # Set column names
    colnames(Eref) <- names_Eref
    colnames(CMD) <- names_CMD
    colnames(ErefR) <- names_ErefR
    Eref <- data.table::as.data.table(Eref)
    CMD <- data.table::as.data.table(CMD)
    ErefR <- data.table::as.data.table(ErefR)
    return(list(Eref = Eref, CMD = CMD, ErefR = ErefR))
  }
  #rd <- calculate_Eref_CMD(xyl_set, Tmax3, Tmin3, Tave3, ppt3) ;head(rd)
  #rd <- calculate_Eref_CMD(xyl_set, Tmax12, Tmin12, Tave12, PPT12) ;head(rd)
  #Eref <- rd$Eref
  # end of Functions =====================================================================<
  
  
  ## pre-set a data.table for var_all ==========================================
  varList_Y=c("MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","bFFP","eFFP","FFP","CMD","CMI","DD_0","DD5","DD_18","DD18","DD1040","EMT","EXT",
              "Eref", "rsds","NFFD", "PAS","RH")
  varList_S=c("Tmax_wt","Tmax_sp","Tmax_sm","Tmax_at","Tmin_wt","Tmin_sp","Tmin_sm","Tmin_at","Tave_wt","Tave_sp","Tave_sm","Tave_at",
              "PPT_wt","PPT_sp","PPT_sm","PPT_at","rsds_wt","rsds_sp","rsds_sm","rsds_at",
              "DD_0_wt","DD_0_sp","DD_0_sm","DD_0_at","DD5_wt","DD5_sp","DD5_sm","DD5_at","DD_18_wt",
              "DD_18_sp","DD_18_sm","DD_18_at","DD18_wt","DD18_sp","DD18_sm","DD18_at","NFFD_wt","NFFD_sp",
              "NFFD_sm","NFFD_at","PAS_wt","PAS_sp","PAS_sm","PAS_at","Eref_wt","Eref_sp","Eref_sm","Eref_at","CMD_wt",
              "CMD_sp","CMD_sm","CMD_at","RH_wt","RH_sp","RH_sm","RH_at","CMI_wt","CMI_sp","CMI_sm","CMI_at")
  varList_PM=c("Tmax01","Tmax02","Tmax03","Tmax04","Tmax05","Tmax06","Tmax07","Tmax08","Tmax09","Tmax10","Tmax11","Tmax12",
               "Tmin01","Tmin02","Tmin03","Tmin04","Tmin05","Tmin06","Tmin07","Tmin08","Tmin09","Tmin10","Tmin11","Tmin12",
               "Tave01","Tave02","Tave03","Tave04","Tave05","Tave06","Tave07","Tave08","Tave09","Tave10","Tave11","Tave12",
               "PPT01","PPT02","PPT03","PPT04","PPT05","PPT06","PPT07","PPT08","PPT09","PPT10","PPT11","PPT12",
               "rsds01","rsds02","rsds03","rsds04","rsds05","rsds06","rsds07","rsds08","rsds09","rsds10","rsds11","rsds12")
  varList_DM=c("DD_0_01","DD_0_02","DD_0_03","DD_0_04","DD_0_05","DD_0_06","DD_0_07","DD_0_08","DD_0_09","DD_0_10", "DD_0_11","DD_0_12",
               "DD5_01","DD5_02","DD5_03","DD5_04","DD5_05","DD5_06","DD5_07","DD5_08","DD5_09","DD5_10","DD5_11","DD5_12",
               "DD_18_01","DD_18_02","DD_18_03","DD_18_04","DD_18_05","DD_18_06","DD_18_07","DD_18_08","DD_18_09","DD_18_10", "DD_18_11","DD_18_12",
               "DD18_01","DD18_02","DD18_03","DD18_04","DD18_05","DD18_06","DD18_07","DD18_08","DD18_09","DD18_10","DD18_11","DD18_12",
               "NFFD01","NFFD02","NFFD03","NFFD04","NFFD05","NFFD06","NFFD07","NFFD08","NFFD09","NFFD10","NFFD11","NFFD12",
               "PAS01","PAS02","PAS03","PAS04","PAS05","PAS06","PAS07","PAS08","PAS09","PAS10","PAS11","PAS12",
               "Eref01","Eref02","Eref03","Eref04","Eref05","Eref06","Eref07","Eref08","Eref09","Eref10","Eref11","Eref12",
               "CMD01","CMD02","CMD03","CMD04","CMD05","CMD06","CMD07","CMD08","CMD09","CMD10","CMD11","CMD12",
               "RH01","RH02","RH03","RH04","RH05","RH06","RH07","RH08","RH09","RH10","RH11","RH12",
               "CMI01","CMI02","CMI03","CMI04","CMI05","CMI06","CMI07","CMI08","CMI09","CMI10","CMI11","CMI12")
  
  var0_pdd = c("PPT01","PPT02","PPT03","PPT04","PPT05","PPT06","PPT07","PPT08","PPT09","PPT10","PPT11","PPT12","PPT_wt","PPT_sp","PPT_sm","PPT_at","MAP","MSP",
               "DD_0_01","DD_0_02","DD_0_03","DD_0_04","DD_0_05","DD_0_06","DD_0_07","DD_0_08","DD_0_09","DD_0_10", "DD_0_11","DD_0_12","DD_0_wt","DD_0_sp","DD_0_sm","DD_0_at","DD_0",
               "DD5_01","DD5_02","DD5_03","DD5_04","DD5_05","DD5_06","DD5_07","DD5_08","DD5_09","DD5_10","DD5_11","DD5_12","DD5_wt","DD5_sp","DD5_sm","DD5_at","DD5",
               "DD61_01","DD61_02","DD61_03","DD61_04","DD61_05","DD61_06","DD61_07","DD61_08","DD61_09","DD61_10","DD61_11","DD61_12")
  var0_ddnp = c("DD_18_01","DD_18_02","DD_18_03","DD_18_04","DD_18_05","DD_18_06","DD_18_07","DD_18_08","DD_18_09","DD_18_10", "DD_18_11","DD_18_12","DD_18_wt","DD_18_sp","DD_18_sm","DD_18_at","DD_18",
                "DD18_01","DD18_02","DD18_03","DD18_04","DD18_05","DD18_06","DD18_07","DD18_08","DD18_09","DD18_10","DD18_11","DD18_12","DD18_wt","DD18_sp","DD18_sm","DD18_at","DD18",
                "NFFD01","NFFD02","NFFD03","NFFD04","NFFD05","NFFD06","NFFD07","NFFD08","NFFD09","NFFD10","NFFD11","NFFD12","NFFD_wt","NFFD_sp","NFFD_sm","NFFD_at","NFFD",
                "PAS01","PAS02","PAS03","PAS04","PAS05","PAS06","PAS07","PAS08","PAS09","PAS10","PAS11","PAS12","PAS_wt","PAS_sp","PAS_sm","PAS_at","PAS")
  var0_ecr = c("Eref01","Eref02","Eref03","Eref04","Eref05","Eref06","Eref07","Eref08","Eref09","Eref10","Eref11","Eref12","Eref_wt","Eref_sp","Eref_sm","Eref_at", "Eref",
               "CMD01","CMD02","CMD03","CMD04","CMD05","CMD06","CMD07","CMD08","CMD09","CMD10","CMD11","CMD12","CMD_wt","CMD_sp","CMD_sm","CMD_at","CMD",
               "RH01","RH02","RH03","RH04","RH05","RH06","RH07","RH08","RH09","RH10","RH11","RH12","RH_wt","RH_sp","RH_sm","RH_at","RH",
               "bFFP","eFFP","FFP","DD1040","ErefR01","ErefR02","ErefR03","ErefR04","ErefR05","ErefR06","ErefR07","ErefR08","ErefR09","ErefR10","ErefR11","ErefR12")
  var1_ttr = c("Tmax01","Tmax02","Tmax03","Tmax04","Tmax05","Tmax06","Tmax07","Tmax08","Tmax09","Tmax10","Tmax11","Tmax12","Tmax_wt","Tmax_sp","Tmax_sm","Tmax_at","MWMT",
               "Tmin01","Tmin02","Tmin03","Tmin04","Tmin05","Tmin06","Tmin07","Tmin08","Tmin09","Tmin10","Tmin11","Tmin12","Tmin_wt","Tmin_sp","Tmin_sm","Tmin_at","MCMT",
               "Tave01","Tave02","Tave03","Tave04","Tave05","Tave06","Tave07","Tave08","Tave09","Tave10","Tave11","Tave12","Tave_wt","Tave_sp","Tave_sm","Tave_at","MAT",
               "rsds01","rsds02","rsds03","rsds04","rsds05","rsds06","rsds07","rsds08","rsds09","rsds10","rsds11","rsds12","rsds_wt","rsds_sp","rsds_sm","rsds_at","rsds",
               "TD","AHM","SHM","EMT","EXT")
  var2_cmi = c("CMI01","CMI02","CMI03","CMI04","CMI05","CMI06","CMI07","CMI08","CMI09","CMI10","CMI11","CMI12","CMI_wt","CMI_sp","CMI_sm","CMI_at","CMI")
  
  var_bsl = c("Tmax01","Tmax02","Tmax03","Tmax04","Tmax05","Tmax06","Tmax07","Tmax08","Tmax09","Tmax10","Tmax11","Tmax12",
              "Tmin01","Tmin02","Tmin03","Tmin04","Tmin05","Tmin06","Tmin07","Tmin08","Tmin09","Tmin10","Tmin11","Tmin12",
              "PPT01","PPT02","PPT03","PPT04","PPT05","PPT06","PPT07","PPT08","PPT09","PPT10","PPT11","PPT12",
              "rsds01","rsds02","rsds03","rsds04","rsds05","rsds06","rsds07","rsds08","rsds09","rsds10","rsds11","rsds12","elev")
  
  length(var0_pdd);length(var0_ddnp);length(var0_ecr);length(var1_ttr);length(var2_cmi);gc()
  DT0_pdd <- data.table(matrix(NA_integer_, nrow=nrow(xyl_input), ncol=length(var0_pdd))); setnames(DT0_pdd, var0_pdd);DT0_pdd;gc()
  DT0_ddnp <- data.table(matrix(NA_integer_, nrow=nrow(xyl_input), ncol=length(var0_ddnp))); setnames(DT0_ddnp, var0_ddnp);DT0_ddnp;gc()
  DT0_ecr <- data.table(matrix(NA_integer_, nrow=nrow(xyl_input), ncol=length(var0_ecr))); setnames(DT0_ecr, var0_ecr);DT0_ecr;gc()
  DT1_ttr <- data.table(matrix(NA_integer_, nrow=nrow(xyl_input), ncol=length(var1_ttr))); setnames(DT1_ttr, var1_ttr);DT1_ttr;gc()
  DT2_cmi <- data.table(matrix(NA_integer_, nrow=nrow(xyl_input), ncol=length(var2_cmi))); setnames(DT2_cmi, var2_cmi);DT2_cmi;gc(); object.size(DT0_pdd)
  DT_bsl  <- data.table(matrix(NA_real_, nrow=nrow(xyl_input), ncol=length(var_bsl))); setnames(DT_bsl, var_bsl);head(DT_bsl);gc(); object.size(DT_bsl)
  rm(var0_pdd, var0_ddnp, var0_ecr,var_bsl);gc();ls()
  
  
  ## baseline data loading and downscaling ========================================================
  #baseline data loading and making stk ---------------
  #terra::terraOptions(memfrac = 0.65)  # Use 50% of available memory
  if(input_valid==T){
    message('Loading baseline data (it takes a minute) ...')
    csv_bs <- data.table(ClimateNAr:::id_bs,
                         load_object(system.file('data/baseline', 'dv_bs.rda', package = 'ClimateNAr'))
    )
    names(csv_bs);dim(csv_bs)
    stk_el  <- terra::rast(csv_bs[, c(1:3), with = FALSE]);stk_el
    stk_tmx <- terra::rast(csv_bs[, c(1:2, 4:15), with = FALSE]);stk_tmx
    stk_tmn <- terra::rast(csv_bs[, c(1:2, 16:27), with = FALSE]);stk_tmn;gc()
    stk_ppt <- terra::rast(csv_bs[, c(1:2, 28:39), with = FALSE]);stk_ppt
    #stk_rad <- terra::rast(csv_bs[, c(1:2, 40:51), with = FALSE]);stk_rad; gc()
    stk_tmxr <- terra::rast(csv_bs[, c(1:2, 52:63), with = FALSE]);stk_tmxr
    stk_tmnr <- terra::rast(csv_bs[, c(1:2, 64:75), with = FALSE]);stk_tmnr; gc()
    stk_pptr <- terra::rast(csv_bs[, c(1:2, 76:87), with = FALSE]);stk_pptr
    #stk_radr <- terra::rast(csv_bs[, c(1:2, 88:99), with = FALSE]);stk_radr
    rm(csv_bs);gc()
  }else{stop('Input file is not valid')}
  
  ## baseline data downscaling  --------------------
  message('Baseline data downscaling ...')
  ## splitting data into sets to process -----------------------------
  setSize = 500000; nSet = floor(nrow(xyl_input)/setSize)+1;nSet
  Set=1
  for(Set in 1:nSet){
    ## extract baseline climate -----------------------------
    row0=(Set-1)*setSize+1;row0
    row2= row0+setSize-1;row2
    if(Set==nSet){row2= nrow(xyl_input);row2} #for the last set
    nrm_el = data.table::as.data.table(terra::extract(stk_el,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_el);dim(nrm_el)
    xyl_set <- data.table::data.table(xyl_input[row0:row2,1:3],nrm_el);names(xyl_set)[1:2]<- c('x','y');head(xyl_set)
    
    nrm_tmx <- data.table::as.data.table(terra::extract(stk_tmx,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_tmx,2);dim(nrm_tmx)
    nrm_tmxr <- data.table::as.data.table(terra::extract(stk_tmxr,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_tmxr,2);dim(nrm_tmxr)
    nrm_tmn  <- data.table::as.data.table(terra::extract(stk_tmn,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_tmn,2);dim(nrm_tmn)
    nrm_tmnr <- data.table::as.data.table(terra::extract(stk_tmnr,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_tmnr,2);dim(nrm_tmnr)
    nrm_ppt  <- data.table::as.data.table(terra::extract(stk_ppt,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_ppt,2);dim(nrm_ppt)
    nrm_pptr  <- data.table::as.data.table(terra::extract(stk_pptr,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_pptr,2);dim(nrm_pptr)
    # nrm_rad  <- data.table::as.data.table(terra::extract(stk_rad,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_rad,2);dim(nrm_rad)
    # nrm_radr  <- data.table::as.data.table(terra::extract(stk_radr,xyl_input[row0:row2,1:2],method="bilinear",ID=F));head(nrm_radr,2);dim(nrm_radr)
    gc();ls()
    
    ## downscaling ----------------
    head(xyl_set);dim(xyl_set)
    el_D <- unlist(xyl_set$elev0 - xyl_set$elev);head(el_D);length(el_D)
    #if(inputFile_class!='Raster'){el_D[is.na(el_D)] <- 0}
    Tmax <- nrm_tmx/100;head(Tmax);dim(Tmax)
    Tmin <- nrm_tmn/100;head(Tmin);dim(Tmin)
    ppt <- data.frame(nrm_ppt);head(ppt);dim(ppt)
    #rad <- nrm_rad/10;head(rad);dim(rad)
    rm(nrm_tmx,nrm_tmn,nrm_ppt);gc()
    
    Tmax_lr <- nrm_tmxr/100000;head(Tmax_lr)
    Tmin_lr <- nrm_tmnr/10000;head(Tmin_lr)
    ppt_lr <- data.frame(nrm_pptr/100000);head(ppt_lr)
    #rad_lr <- nrm_radr/100000;head(rad_lr)
    rm(nrm_tmxr,nrm_tmnr,nrm_pptr);gc()
    
    # Elevation adjustment------------------------------------------
    Tmax2 <- Tmax + Tmax_lr*el_D; names(Tmax2) <- gsub('t','T',names(Tmax2));head(Tmax2)
    Tmin2 <- Tmin + Tmin_lr*el_D;names(Tmin2) <- gsub('t','T',names(Tmin2));head(Tmin2)
    Tave2 <- (Tmax2+Tmin2)/2; names(Tave2) <- paste0('Tave',substr(101:112,2,3));head(Tave2)
    ppt2  <- ppt + ppt_lr *el_D; names(ppt2) <- gsub('prec','PPT',names(ppt2));head(ppt2)
    ppt2[ppt2<0] <- 0
    #rad2 <- rad + rad_lr*el_D;head(rad2)
    #rad2[rad2 < -100] <- NA
    rm(Tmax,Tmax_lr,Tmin,Tmin_lr,ppt,ppt_lr,el_D);gc(); ls()
    
    # rsds --------------
    ecd <- calculate_Eref_CMD(xyl_set, Tmax2, Tmin2, Tave2, ppt2) ;head(ecd)
    ErefR <- ecd$ErefR;head(ErefR)
    rh0 <- calculate_RH(Tmax2, Tmin2);head(rh0)
    td03 = (Tmax2$Tmax03 - Tmin2$Tmin03);head(td03)
    rsds <- data.table(rsds01 = 29.99 + (-0.716 * xyl_set$y) + (-0.185 * ErefR$ErefR01) + (0.00326 * xyl_set$elev) + (0.00406 * xyl_set$y^2) + (0.000787 * ErefR$ErefR01^2) + (-0.000000348 * xyl_set$elev^2) + (0.00357 * xyl_set$y * ErefR$ErefR01) + (-0.0000480 * xyl_set$y * xyl_set$elev) + (-0.0000278 * ErefR$ErefR01 * xyl_set$elev) + (0.00000150 * xyl_set$y * ErefR$ErefR01 * xyl_set$elev))
    rsds$rsds02 = -8.66 + (0.206 * xyl_set$y) + (-0.491 * xyl_set$x) + (0.473 * rh0$RH02) + (0.000856 * xyl_set$y^2) + (0.000169 * xyl_set$x^2) + (-0.000741 * rh0$RH02^2) + (0.00893 * xyl_set$y * xyl_set$x) + (-0.00718 * xyl_set$y * rh0$RH02) + (0.00703 * xyl_set$x * rh0$RH02) + (-0.000119 * xyl_set$y * xyl_set$x * rh0$RH02)
    rsds$rsds03 = 53.4 + (-0.613 * xyl_set$y) + (-3.72 * td03) + (0.432 * xyl_set$x) + (-0.00168 * xyl_set$y^2) + (-0.0112 * td03^2) + (0.000203 * xyl_set$x^2) + (0.0680 * xyl_set$y * td03) + (-0.00620 * xyl_set$y * xyl_set$x) + (-0.0472 * td03 * xyl_set$x) + (0.000769 * xyl_set$y * td03 * xyl_set$x)
    rsds$rsds04 = -57.2 + (1.18 * xyl_set$y) + (1.18 * rh0$RH04) + (-1.00 * xyl_set$x) + (-0.000156 * xyl_set$y^2) + (-0.000715 * rh0$RH04^2) + (0.000147 * xyl_set$x^2) + (-0.0179 * xyl_set$y * rh0$RH04) + (0.0160 * xyl_set$y * xyl_set$x) + (0.0142 * rh0$RH04 * xyl_set$x) + (-0.000220 * xyl_set$y * rh0$RH04 * xyl_set$x)
    rsds$rsds05 = 21.9 + (-0.0742 * ErefR$ErefR05) + (0.00554 * xyl_set$elev) + (-0.00834 * xyl_set$x) + (0.000243 * ErefR$ErefR05^2) + (0.000000241 * xyl_set$elev^2) + (-0.000211 * xyl_set$x^2) + (-0.000138 * ErefR$ErefR05 * xyl_set$elev) + (-0.000472 * ErefR$ErefR05 * xyl_set$x) + (0.0000479 * xyl_set$elev * xyl_set$x) + (-0.00000135 * ErefR$ErefR05 * xyl_set$elev * xyl_set$x)
    rsds$rsds06 = 83.5 + (-1.32 * rh0$RH06) + (-0.0587 * xyl_set$x) + (-1.48 * xyl_set$y) + (0.00314 * rh0$RH06^2) + (-0.000728 * xyl_set$x^2) + (-0.00205 * xyl_set$y^2) + (-0.00183 * rh0$RH06 * xyl_set$x) + (0.0240 * rh0$RH06 * xyl_set$y) + (-0.00745 * xyl_set$x * xyl_set$y) + (0.000110 * rh0$RH06 * xyl_set$x * xyl_set$y)
    rsds$rsds07 = 7.70 + (-0.00986 * ErefR$ErefR07) + (-0.198 * xyl_set$x) + (0.0158 * xyl_set$elev) + (0.000164 * ErefR$ErefR07^2) + (-0.000995 * xyl_set$x^2) + (0.000000309 * xyl_set$elev^2) + (-0.000159 * ErefR$ErefR07 * xyl_set$x) + (-0.000190 * ErefR$ErefR07 * xyl_set$elev) + (0.000143 * xyl_set$x * xyl_set$elev) + (-0.00000177 * ErefR$ErefR07 * xyl_set$x * xyl_set$elev)
    rsds$rsds08 = 5.79 + (0.0532 * ErefR$ErefR08) + (-0.0923 * xyl_set$x) + (0.00843 * xyl_set$elev) + (-0.000157 * ErefR$ErefR08^2) + (-0.000654 * xyl_set$x^2) + (0.000000281 * xyl_set$elev^2) + (-0.000533 * ErefR$ErefR08 * xyl_set$x) + (-0.000138 * ErefR$ErefR08 * xyl_set$elev) + (0.0000691 * xyl_set$x * xyl_set$elev) + (-0.00000126 * ErefR$ErefR08 * xyl_set$x * xyl_set$elev)
    rsds$rsds09 = 4.90 + (0.102 * ErefR$ErefR09) + (-0.000119 * xyl_set$elev) + (-0.00471 * xyl_set$x) + (-0.000457 * ErefR$ErefR09^2) + (0.000000371 * xyl_set$elev^2) + (-0.000133 * xyl_set$x^2) + (-0.0000745 * ErefR$ErefR09 * xyl_set$elev) + (-0.000686 * ErefR$ErefR09 * xyl_set$x) + (-0.00000248 * xyl_set$elev * xyl_set$x) + (-0.000000721 * ErefR$ErefR09 * xyl_set$elev * xyl_set$x)
    rsds$rsds10 = -0.191 + (0.0901 * ErefR$ErefR10) + (0.216 * xyl_set$y) + (0.000974 * xyl_set$elev) + (-0.00262 * xyl_set$y^2) + (0.000000234 * xyl_set$elev^2) + (0.000349 * ErefR$ErefR10 * xyl_set$y) + (-0.0000175 * ErefR$ErefR10 * xyl_set$elev) + (-0.0000160 * xyl_set$y * xyl_set$elev) + (0.000000964 * ErefR$ErefR10 * xyl_set$y * xyl_set$elev)
    rsds$rsds11= 11.6 + (-0.184 * xyl_set$y) + (-0.0492 * ErefR$ErefR11) + (0.00265 * xyl_set$elev) + (0.000313 * xyl_set$y^2) + (0.000468 * ErefR$ErefR11^2) + (-0.000000115 * xyl_set$elev^2) + (0.00272 * xyl_set$y * ErefR$ErefR11) + (-0.0000391 * xyl_set$y * xyl_set$elev) + (-0.0000282 * ErefR$ErefR11 * xyl_set$elev) + (0.00000144 * xyl_set$y * ErefR$ErefR11 * xyl_set$elev)
    rsds$rsds12 = 20.5 + (-0.482 * xyl_set$y) + (-0.110 * ErefR$ErefR12) + (0.00364 * xyl_set$elev) + (0.00259 * xyl_set$y^2) + (0.000584 * ErefR$ErefR12^2) + (-0.000000274 * xyl_set$elev^2) + (0.00316 * xyl_set$y * ErefR$ErefR12) + (-0.0000542 * xyl_set$y * xyl_set$elev) + (-0.0000293 * ErefR$ErefR12 * xyl_set$elev) + (0.00000130 * xyl_set$y * ErefR$ErefR12 * xyl_set$elev)
    rsds[rsds<0] <- 0
    rsds <- round(rsds,1);   head(rsds)
    
    # add toDT_bsl ------------------
    DT_bsl[row0:row2,c(paste0('Tmax0',1:9),'Tmax10','Tmax11','Tmax12') := round(Tmax2,2)];DT_bsl
    DT_bsl[row0:row2,c(paste0('Tmin0',1:9),'Tmin10','Tmin11','Tmin12') := round(Tmin2,2)];DT_bsl
    DT_bsl[row0:row2,c(paste0('PPT0',1:9),'PPT10','PPT11','PPT12') := round(ppt2,0)];DT_bsl
    #DT_bsl[row0:row2,c(paste0('Rad0',1:9),'Rad10','Rad11','Rad12') := round(rad2,1)];DT_bsl
    DT_bsl[row0:row2,names(rsds) := rsds];DT_bsl
    DT_bsl[row0:row2,'elev' := round(nrm_el,0)];DT_bsl
    rm(Tmax2,Tmin2,Tave2,ppt2,nrm_el,ecd,ErefR,rh0,td03,rsds);gc()
  }
  rm(stk_el,stk_tmx,stk_tmxr,stk_tmn,stk_tmnr,stk_ppt,stk_pptr);gc()
  
  ## Reference, historical and/or GCM ===============================================================>
  ## Reference, historical and/or GCM ===============================================================>
  message('Generating data for periods ....')
  # periodList= c('Normal_1961_1990.nrm','Year_1902.ann','8GCMs_ensemble_ssp245_2041-2070.gcm')
  # periodList= c('8GCMs_ensemble_ssp245_2041-2070.gcm')
  # periodList= c('Normal_1961_1990.nrm')
  # periodList= c('Normal_1971_2000.nrm')
  # periodList= c('Normal_1971_1980.dcd')
  # periodList= 'Year_1902.ann'
  # periodList= 1901:1905
  # periodList= '8GCM_ssp245_2301'
  # periodList= '8GCM_ssp245_2031:2040'
  
  GCM_ts = FALSE;yr=0
  periodList2 <- periodList;periodList2
  if(nchar(periodList[1])==21){ #gcm time series
    ptext = (substr(periodList,13,21));ptext
    periodList2 <- eval(parse(text = ptext));periodList2
    GCM_ts = TRUE
  };periodList2
  
  prd = periodList2[1];prd
  for(prd in periodList2){
    # determine period_class: 'His', 'His_Ann', 'GCM', 'GCM_Ann' or 'Ref' -------------
    prd_class='Ref' #default
    
    #historical ---------------
    if(grepl('.nrm|.dcd',prd)){ # historical period
      prd_class='His_prd'
    };prd_class
    if(grepl('Year_',prd)){ #individual historical years
      prd_class='His_Ann'
      yr = as.numeric(gsub("[^0-9]","",prd))
    };prd_class
    if(is.numeric(prd)&GCM_ts==FALSE){ # his_ts
      prd_class='His_Ann'
      yr = prd
      prd = paste0('Year_',prd,'.ann');prd
    };prd_class
    
    #GCMs ----------------
    if(grepl('.gcm',prd)){prd_class='GCM'};prd_class #gcm normals
    if(grepl('8GCM_ssp',prd)){ #gcm individual years
      prd_class='GCM_Ann'
      prd = paste0(substr(prd,1,11),'/', as.numeric(substr(prd,13,16)),'.gcm');prd
    };prd_class
    if(is.numeric(prd)&GCM_ts==TRUE){ #gcm_ts
      prd_class='GCM_Ann'
      prd = paste0(substr(periodList,1,11),'/',prd,'.gcm');prd
    };prd_class;prd
    
    #referecne period -------
    if(prd =='Normal_1961_1990.nrm'){prd_class='Ref'};prd_class
    
    
    ## splitting data into sets to process -----------------------------
    stk_h_open =0; prd_prev = 0; stk_gcm_open=0 #dataset openning status
    setSize = 500000; nSet = floor(nrow(xyl_input)/setSize)+1;nSet
    Set=1
    for(Set in 1:nSet){
      ## extract baseline climate -----------------------------
      row0=(Set-1)*setSize+1;row0
      row2= row0+setSize-1;row2
      if(Set==nSet){row2= nrow(xyl_input);row2} #for the last set
      elev <- DT_bsl$elev[row0:row2];head(elev)
      xyl_set <- data.table(xyl_input[row0:row2,1:3],elev);names(xyl_set)[1:2]<- c('x','y');head(xyl_set)
      Tmax2 <- DT_bsl[row0:row2,1:12];head(Tmax2)
      Tmin2 <- DT_bsl[row0:row2,13:24];head(Tmin2)
      ppt2 <-  DT_bsl[row0:row2,25:36];head(ppt2)
      #rad2 <-  DT_bsl[row0:row2,37:48];head(rad2)
      rsds2 <- DT_bsl[row0:row2,37:48];head(rsds2)
      
      ## for reference period -------------------------
      if(prd_class=='Ref'){
        Tmax3 = round(Tmax2,1)
        Tmin3 = round(Tmin2,1)
        ppt3 = round(ppt2,0)
        #rad3 = round(rad2,1)
        rsds3 <- rsds2
        Tave3 <- round((Tmax3+Tmin3)/2, 1); names(Tave3) <- paste0('Tave',substr((101:112),2,3));head(Tave3)
        Tmax12 = Tmax3[,12] #for seasonal vars
        Tmin12 = Tmin3[,12]
        Tave12 = Tave3[,12]
        PPT12 = ppt3[,12]
      }#end of Ref
      
      # historical stk_hist -------
      if(grepl('His',prd_class)){
        #load historical data ---
        if(stk_h_open !=1|prd != prd_prev){
          hFile = system.file('data/historical', paste0(prd,'.rda'),package="ClimateNAr");hFile
          csv_h <- load_object(hFile);head(csv_h)
          id_h <- load_object(system.file('data/historical','na_cru_id.rda',package='ClimateNAr'));head(id_h)
          xy_h = id_h[,1:2];head(xy_h)
          xyv_h <- data.table(xy_h, csv_h);head(xyv_h)
          stk_h <- terra::rast(xyv_h);stk_h
          stk_h_open =1
          prd_prev = prd
          
          #historical previous year
          if(prd_class=='His_Ann'& yr>1901){
            hFile_py = system.file('data/historical', paste0('Year_',yr-1,'.ann.rda'),package="ClimateNAr");hFile_py
            csv_h_py <- load_object(hFile_py);head(csv_h_py)
            xyv_h_py <- data.table(xy_h, csv_h_py);head(xyv_h_py)
            stk_h_py <- terra::rast(xyv_h_py);stk_h_py
          }
          
          #his normal
          if(prd_class=='His_Ann'){
            if(yr < 1921){nrmF0 = "Normal_1901_1930.nrm"}
            if(yr > 1920 & yr < 1931){nrmF0 = "Normal_1911_1940.nrm"}
            if(yr > 1930 & yr < 1941){nrmF0 = "Normal_1921_1950.nrm"}
            if(yr > 1940 & yr < 1951){nrmF0 = "Normal_1931_1960.nrm"}
            if(yr > 1950 & yr < 1961){nrmF0 = "Normal_1941_1970.nrm"}
            if(yr > 1960 & yr < 1971){nrmF0 = "Normal_1951_1980.nrm"}
            if(yr > 1970 & yr < 1981){nrmF0 = "Normal_1961_1990.nrm"}
            if(yr > 1980 & yr < 1991){nrmF0 = "Normal_1971_2000.nrm"}
            if(yr > 1990 & yr < 2001){nrmF0 = "Normal_1981_2010.nrm"}
            if(yr > 2000){ nrmF0 = "Normal_1991_2020.nrm"}
            nrmF = system.file('data/historical', paste0(nrmF0,'.rda'),package="ClimateNAr");nrmF
            cruNrm <- load_object(nrmF);head(cruNrm)
            xyv_cruNrm <- data.table(xy_h,cruNrm)
            cruNrm_stk <- terra::rast(xyv_cruNrm);cruNrm_stk
          }
        }
        
        # process hist data ---
        dat_h <- data.table(terra::extract(stk_h,xyl_set[,1:2],method="bilinear"));head(dat_h)
        Tmax_h <- dat_h[,14:25];head(Tmax_h)
        Tmin_h <- dat_h[,2:13];head(Tmin_h)
        ppt_h <- dat_h[,26:37];head(ppt_h)
        #rad_h <- dat_h[,38:49];head(rad_h)
        rm(dat_h);gc()
        
        # adding hist data ----------------------
        Tmax3 <- round(Tmax2 + Tmax_h,1); head(Tmax3)
        Tmin3 <- round(Tmin2 + Tmin_h,1); head(Tmin3)
        ppt3 <- round(ppt2 *(100+ppt_h)/100,0); head(ppt3)
        ppt3[ppt3<0] <- 0
        #rad3 <- round(rad2 + rad_h,1); head(rad3)
        #rad3[rad3 < -1000] <- NA;head(rad3)
        Tave3 = round((Tmax3+Tmin3)/2,1);  names(Tave3) <- paste0('Tave',substr(101:112,2,3));head(Tave3)
        Tmax12 = Tmax3[,12] #for seasonal vars without conidering previous year
        Tmin12 = Tmin3[,12]
        Tave12 = Tave3[,12]
        PPT12 = ppt3[,12]
        rsds3 = rsds2
        
        # previous year
        if(prd_class=='His_Ann'& yr>1901){
          dat_h_py <- data.table(terra::extract(stk_h_py,xyl_set[,1:2],method="bilinear"));head(dat_h_py)
          Tmax_h_py <- dat_h_py[,14:25];head(Tmax_h_py)
          Tmin_h_py <- dat_h_py[,2:13];head(Tmin_h_py)
          ppt_h_py <- dat_h_py[,26:37];head(ppt_h_py)
          # rad_h_py <- dat_h_py[,38:49];head(rad_h_py)
          rm(dat_h_py);gc()
          
          # adding hist data ----------------------
          Tmax3_py <- Tmax2 + Tmax_h_py; head(Tmax3_py)
          Tmin3_py <- Tmin2 + Tmin_h_py; head(Tmin3_py)
          ppt3_py <- ppt2 *(100+ppt_h_py)/100; head(ppt3_py)
          ppt3_py[ppt3_py<0] <- 0
          #rad3_py <- rad2 + rad_h_py; head(rad3_py)
          #rad3_py[rad3_py < -5000] <- NA;head(rad3_py)
          Tave3_py = (Tmax3_py+Tmin3_py)/2;names(Tave3_py) <- paste0('Tave',substr(101:112,2,3));head(Tave3)
          Tmax12 = Tmax3_py[,12] #for seasonal vars
          Tmin12 = Tmin3_py[,12]
          Tave12 = Tave3_py[,12];head(Tave12)
          PPT12 = ppt3_py[,12]
        }
        
        #nrm for emt, ext
        if(prd_class=='His_Ann'){
          cruNrm_ex <- data.table(terra::extract(cruNrm_stk,xyl_set[,1:2],method="bilinear"));head(cruNrm_ex)
          Tmax_cruNrm <- cruNrm_ex[,14:25];head(Tmax_cruNrm)
          Tmin_cruNrm <- cruNrm_ex[,2:13];head(Tmin_cruNrm)
          Tmax3_cruNrm <- Tmax2 + Tmax_cruNrm; head(Tmax3_cruNrm)
          Tmin3_cruNrm <- Tmin2 + Tmin_cruNrm; head(Tmin3_cruNrm)
        }
      }#end of stk_h
      
      
      # GCM stk_gcm  --------------------------------------------------
      if(prd_class=='GCM'|prd_class=='GCM_Ann'){
        if(stk_gcm_open !=1|prd != prd_prev){
          if(prd_class=='GCM'){gcmFile = system.file('data/GCMs',paste0(prd,'.rda'), package='ClimateNAr');gcmFile}
          if(prd_class=='GCM_Ann'){gcmFile = system.file('data/GCMs/Annual',paste0(prd,'.rda'), package='ClimateNAr');gcmFile}
          csv_gcm <- load_object(gcmFile);head(csv_gcm)
          id_gcm = load_object(system.file('data/GCMs','na_gcm_id.rda',package='ClimateNAr'));head(csv_gcm)
          xy_gcm = id_gcm[,1:2];head(xy_gcm)
          xyv_gcm <- data.table(xy_gcm, csv_gcm);head(xyv_gcm)
          stk_gcm <- terra::rast(xyv_gcm);stk_gcm
          stk_gcm_open =1
          prd_prev = prd
        }
        
        dat_gcm <- data.table(terra::extract(stk_gcm,xyl_set[,1:2],method="bilinear"));head(dat_gcm);names(dat_gcm)
        Tmax_gcm <- dat_gcm[,2:13];head(Tmax_gcm)
        Tmin_gcm <- dat_gcm[,14:25];head(Tmin_gcm)
        ppt_gcm <- dat_gcm[,26:37];head(ppt_gcm)
        rad_gcm <- dat_gcm[,38:49];head(rad_gcm)
        
        # adding GCM data --------------------------
        Tmax3 <- round(Tmax2 + Tmax_gcm, 1); head(Tmax3)
        Tmin3 <- round(Tmin2 + Tmin_gcm,1); head(Tmin3)
        ppt3 <- round(ppt2 *(100+ppt_gcm)/100,0); head(ppt3)
        ppt3[ppt3<0] <- 0
        #rad3 <- round(rad2 *(100+rad_gcm)/100,1); head(rad3)
        rsds3 <- round(rsds2 *(100+rad_gcm)/100,1); head(rsds3)
        Tave3 = round((Tmax3+Tmin3)/2,1);  names(Tave3) <- paste0('Tave',substr(101:112,2,3));head(Tave3)
        Tmax12 = Tmax3[,12] #for seasonal vars
        Tmin12 = Tmin3[,12]
        Tave12 = Tave3[,12]
        PPT12 = ppt3[,12]
      }#end of stk_gcm
      rm(Tmax2,Tmin2,ppt2,rsds2);gc();ls()
      ## end of ref, hist or GCM ===================-------------------------------<
      
      
      ##end of period specific ============================================
      ## Primary monthly variables collection by DT =================
      #rad3[rad3<0] <- 0
      ppt3[ppt3<0] <- 0
      rsds3[rsds3<0] <- 0
      rsds3$rsds_wt <- (rsds3[,1]+rsds3[,1]+rsds3[,12])/3
      rsds3$rsds_sp <- rowMeans(rsds3[,3:5])
      rsds3$rsds_sm <- rowMeans(rsds3[,6:8])
      rsds3$rsds_at <- rowMeans(rsds3[,9:11])
      rsds3$rsds <- rowMeans(rsds3[,1:12])
      rsds3 <- round(rsds3,1)
      
      DT1_ttr[row0:row2,c(paste0('Tmax0',1:9),'Tmax10','Tmax11','Tmax12') := Tmax3*10];DT1_ttr
      DT1_ttr[row0:row2,c(paste0('Tmin0',1:9),'Tmin10','Tmin11','Tmin12') := Tmin3*10];DT1_ttr
      DT1_ttr[row0:row2,names(Tave3) := Tave3*10];DT1_ttr
      #DT1_ttr[row0:row2,c(paste0('Rad0',1:9),'Rad10','Rad11','Rad12')  := rad3*10];DT1_ttr
      DT1_ttr[row0:row2,names(rsds3)  := rsds3*10];DT1_ttr
      DT0_pdd[row0:row2,c(paste0('PPT0',1:9),'PPT10','PPT11','PPT12') := ppt3];DT0_pdd
      
      
      ##calculated variables (to all)=====================================
      #annual
      yVar <- data.table(MAT=round(rowMeans(Tave3),1));head(yVar)
      yVar$MAP <- round(rowSums(ppt3),0);head(yVar)
      yVar$MWMT <- round(apply(Tave3, 1, max),1);head(yVar)
      yVar$MCMT <- round(apply(Tave3, 1, min),1);head(yVar)
      yVar$TD <- yVar$MWMT - yVar$MCMT;head(yVar)
      yVar$MSP <- round(rowSums(ppt3[,5:9]),0);head(yVar)
      yVar$AHM <- round((yVar$MAT + 10) / (yVar$MAP / 1000),1); head(yVar)
      yVar$SHM <- round(yVar$MWMT / (yVar$MSP / 1000),1);head(yVar)
      # yVar$MAR <- round(rowMeans(rad3),0);head(yVar);object.size(yVar)
      DT0_pdd[row0:row2,names(yVar)[c(2,6)] := yVar[,c(2,6)]];DT0_pdd
      DT1_ttr[row0:row2,names(yVar)[c(1,3:5,7:8)] := round(yVar[,c(1,3:5,7:8)]*10,0)];DT1_ttr
      gc() #yVar is used later
      
      
      #seasonal
      sVar <- data.table(Tmax_wt = round(rowMeans(Tmax3[,c(1,2,12)]),1));head(sVar)
      sVar$Tmax_sp <- round(rowMeans(Tmax3[,c(3:5)]),1);head(sVar)
      sVar$Tmax_sm <- round(rowMeans(Tmax3[,c(6:8)]),1);head(sVar)
      sVar$Tmax_at <- round(rowMeans(Tmax3[,c(9:11)]),1);head(sVar)
      sVar$Tmin_wt <- round(rowMeans(Tmin3[,c(1,2,12)]),1);head(sVar)
      sVar$Tmin_sp <- round(rowMeans(Tmin3[,c(3:5)]),1);head(sVar)
      sVar$Tmin_sm <- round(rowMeans(Tmin3[,c(6:8)]),1);head(sVar)
      sVar$Tmin_at <- round(rowMeans(Tmin3[,c(9:11)]),1);head(sVar)
      sVar$Tave_wt <- round(rowMeans(Tave3[,c(1,2,12)]),1);head(sVar)
      sVar$Tave_sp <- round(rowMeans(Tave3[,c(3:5)]),1);head(sVar)
      sVar$Tave_sm <- round(rowMeans(Tave3[,c(6:8)]),1);head(sVar)
      sVar$Tave_at <- round(rowMeans(Tave3[,c(9:11)]),1);head(sVar)
      sVar$PPT_wt <- round(rowSums(ppt3[,c(1,2,12)]),0);head(sVar)
      sVar$PPT_sp <- round(rowSums(ppt3[,c(3:5)]),0);head(sVar)
      sVar$PPT_sm <- round(rowSums(ppt3[,c(6:8)]),0);head(sVar)
      sVar$PPT_at <- round(rowSums(ppt3[,c(9:11)]),0);head(sVar)
      #sVar$Rad_wt <- round(rowMeans(rad3[,c(1,2,12)]),1);head(sVar)
      #sVar$Rad_sp <- round(rowMeans(rad3[,c(3:5)]),1);head(sVar)
      #sVar$Rad_sm <- round(rowMeans(rad3[,c(6:8)]),1);head(sVar)
      #sVar$Rad_at <- round(rowMeans(rad3[,c(9:11)]),1);head(sVar)
      
      #consider pre-year
      if(prd_class=='His_Ann'& yr>1901){
        sVar$Tmax_wt <- round((Tmax3[,1]+Tmax3[,2]+Tmax12)/3,1);head(sVar$Tmax_wt)
        sVar$Tmin_wt <- round((Tmin3[,1]+Tmin3[,2]+Tmin12)/3,1);head(sVar$Tmin_wt)
        sVar$Tave_wt = round((sVar$Tmax_wt+sVar$Tmin_wt)/2,1);head(sVar$Tave_wt)
        sVar$PPT_wt <- round((ppt3[,1]+ppt3[,2]+PPT12),0);head(sVar$PPT_wt)
        # sVar$Rad_wt <- round((rad3[,1]+rad3[,2]+rad3_py[,12])/3,0);head(sVar$Rad_wt)
      }
      
      names(sVar)
      DT0_pdd[row0:row2,names(sVar)[c(13:16)] := sVar[,c(13:16)]];DT0_pdd
      DT1_ttr[row0:row2,names(sVar)[c(1:12)] := sVar[,c(1:12)]*10];DT1_ttr
      rm(sVar);gc()
      
      
      #Derived climate variables ------------------>
      #DD_0 --------
      dd_0 <- data.table(DD_0_01=ifelse(Tave3[[1]] >= -9,381.2279 / (1 + exp(-(Tave3[[1]] - -5.46) / -3.5879)),-30.7664 * Tave3[[1]] + 5));head(dd_0)
      dd_0$DD_0_02=ifelse(Tave3[[2]] >= -10,357.1408/ (1 + exp(-(Tave3[[2]] - -5.69) / -3.6034)),-27.7799 * Tave3[[2]] + 5);head(dd_0)
      dd_0$DD_0_03=ifelse(Tave3[[3]] >= -9.5,358.7336/ (1 + exp(-(Tave3[[3]] - -4.97) / -3.3081)),-30.8261 * Tave3[[3]] + 5);head(dd_0)
      dd_0$DD_0_04=ifelse(Tave3[[4]] >= -7, 274.6697/ (1 + exp(-(Tave3[[4]] - -3.59) / -2.6814)),-29.5144 * Tave3[[4]] + 10);head(dd_0)
      dd_0$DD_0_05=ifelse(Tave3[[5]] >= -5, 169.7686/ (1 + exp(-(Tave3[[5]] - -1.23) / -1.8477)),-29.0289 * Tave3[[5]] + 20);head(dd_0)
      dd_0$DD_0_06= 167.0972/ (1 + exp(-(Tave3[[6]] - -2.12) / -1.9556));head(dd_0)
      dd_0$DD_0_07= 134.2072/ (1 + exp(-(Tave3[[7]] - -1.54) / -1.7562));head(dd_0)
      dd_0$DD_0_08= 101.3172/ (1 + exp(-(Tave3[[8]] - -0.96) / -1.5569));head(dd_0)
      dd_0$DD_0_09=ifelse(Tave3[[9]] >= -4.5, 262.5633/ (1 + exp(-(Tave3[[9]] - -3.91) / -2.3346)),-32.6652 * Tave3[[9]] + 5);head(dd_0)
      dd_0$DD_0_10=ifelse(Tave3[[10]] >= -6.5, 246.0955/ (1 + exp(-(Tave3[[10]] - -2.86) / -2.4705)),-30.2925 * Tave3[[10]] + 5);head(dd_0)
      dd_0$DD_0_11=ifelse(Tave3[[11]] >= -7.1, 280.4533/ (1 + exp(-(Tave3[[11]] - -3.56) / -2.8379)),-29.5002 * Tave3[[11]] + 5);head(dd_0)
      dd_0$DD_0_12=ifelse(Tave3[[12]] >= -10, 403.6482/ (1 + exp(-(Tave3[[12]] - -5.95) / -3.6009)),-30.7232 * Tave3[[12]] + 5);head(dd_0)
      dd_0$DD_0_py12=ifelse(Tave12 >= -10, 403.6482/ (1 + exp(-(Tave12[[1]] - -5.95) / -3.6009)),-30.7232 * Tave12[[1]] + 5);head(dd_0)
      
      dd_0$DD_0_wt <- rowSums(dd_0[,c(1,2,13)]);head(dd_0)
      dd_0$DD_0_sp <- rowSums(dd_0[,c(3:5)]);head(dd_0)
      dd_0$DD_0_sm <- rowSums(dd_0[,c(6:8)]);head(dd_0)
      dd_0$DD_0_at <- rowSums(dd_0[,c(9:11)]);head(dd_0)
      dd_0$DD_0 <- rowSums(dd_0[,c(1:12)]);head(dd_0)
      dd_0 <- round(dd_0);head(dd_0)
      
      DT0_pdd[row0:row2,names(dd_0)[c(1:12,14:18)] := dd_0[,c(1:12,14:18)]];DT0_pdd
      rm(dd_0);gc()
      
      #DD5 -------------------------------------------------
      dd5 <- data.table(DD5_01=ifelse(Tave3[[1]] <=	12,	337.5699	/(1+exp(-(Tave3[[1	]]-	10.02	)/	3.3586)),	30.0966	*Tave3[[	1	]]+	-140));head(dd5)
      dd5$DD5_02=ifelse(Tave3[[2]] <=	12,	302.8633/(1+exp(-(Tave3[[2]]-	10.01	)/	3.3086)),	28.0429	*Tave3[[2]]+	-140);head(dd5)
      dd5$DD5_03=ifelse(Tave3[[3]] <=	12,	363.8180/(1+exp(-(Tave3[[3]]-	10.50	)/	3.4363)),	30.1222	*Tave3[[3]]+	-140);head(dd5)
      dd5$DD5_04=ifelse(Tave3[[4]] <=	12,	339.8059/(1+exp(-(Tave3[[4]]-	10.35	)/	3.2628)),	29.4187	*Tave3[[4]]+	-140);head(dd5)
      dd5$DD5_05=ifelse(Tave3[[5]] <=	12,	327.1587/(1+exp(-(Tave3[[5]]-	10.05	)/	2.7967)),	30.1473	*Tave3[[5]]+	-140);head(dd5)
      dd5$DD5_06=ifelse(Tave3[[6]] <=	13,	370.5585/(1+exp(-(Tave3[[6]]-	11.13	)/	3.1320)),	29.9647	*Tave3[[6]]+	-150);head(dd5)
      dd5$DD5_07=ifelse(Tave3[[7]]<=	15,	410.0218/(1+exp(-(Tave3[[7]]-	11.63	)/	3.1190)),	30.7456	*Tave3[[7]]+	-150);head(dd5)
      dd5$DD5_08=ifelse(Tave3[[8]]<=	15,	412.2794/(1+exp(-(Tave3[[8]]-	11.66	)/	3.1271)),	30.7429	*Tave3[[8]]+	-150);head(dd5)
      dd5$DD5_09=ifelse(Tave3[[9]]<=	13,	342.8546/(1+exp(-(Tave3[[9]]-	10.51	)/	2.9657)),	29.7243	*Tave3[[9]]+	-145);head(dd5)
      dd5$DD5_10=ifelse(Tave3[[10]]<=	12,	344.9987/(1+exp(-(Tave3[[10]]-10.26	)/	3.1937)),	30.4110	*Tave3[[10]]+	-145);head(dd5)
      dd5$DD5_11=ifelse(Tave3[[11]]<=	11,	304.7169/(1+exp(-(Tave3[[11]]-9.59	)/	3.1492)),	29.4263	*Tave3[[	11]]+	-140);head(dd5)
      dd5$DD5_12=ifelse(Tave3[[12]]<=	12,	341.0866/(1+exp(-(Tave3[[12]]-10.09	)/	3.2945)),	30.1269	*Tave3[[	12]]+	-140);head(dd5)
      dd5$dd5_py12=ifelse(Tave12[[1]] <=	12,	341.0866/(1+exp(-(Tave12[[1]]-10.09	)/	3.2945)),	30.1269	*Tave12[[1]]+	-140);head(dd5)
      
      ##spatial adjustments ---
      dd5$DD5_01= dd5$DD5_01 + 80.85944122 + -4.292331735 *xyl_set[[2]] + 0.02492521 *xyl_set[[2]]^2 + -0.012754117 *xyl_set[[1]] +  -0.004799503 *xyl_set[[1]]^2 + 1.693255515 *yVar$TD+ -0.026179812 *yVar$TD^2 + -0.017980072 *xyl_set[[2]]*xyl_set[[1]];head(dd5)
      dd5$DD5_02= dd5$DD5_02 + 80.50326623 + -3.375663498 *xyl_set[[2]] + 0.018483736 *xyl_set[[2]]^2 + 0.322029786 *xyl_set[[1]] +  -0.001649449 *xyl_set[[1]]^2 + 1.42412211 *yVar$TD+ -0.01633481 *yVar$TD^2 + -0.012036406 *xyl_set[[2]]*xyl_set[[1]]
      dd5$DD5_03= dd5$DD5_03 + 62.11479743 + -2.544541399 *xyl_set[[2]] + 0.006807518 *xyl_set[[2]]^2 + 0.239920188 *xyl_set[[1]] +  -0.001773369 *xyl_set[[1]]^2 + 1.014172646 *yVar$TD+ 0.002752136 *yVar$TD^2 + -0.010876854 *xyl_set[[2]]*xyl_set[[1]]
      dd5$DD5_04= dd5$DD5_04 + -15.81598921 + -0.690307068 *xyl_set[[2]] + 2.85223E-05 *xyl_set[[2]]^2 + -0.739577834 *xyl_set[[1]] +  -0.004438255 *xyl_set[[1]]^2 + -0.474213942 *yVar$TD+ 0.026258239 *yVar$TD^2 + -0.00240618 *xyl_set[[2]]*xyl_set[[1]]
      dd5$DD5_11= dd5$DD5_11 + 23.28281098 + -1.699943102 *xyl_set[[2]] + 0.00415848 *xyl_set[[2]]^2 + -0.129591771 *xyl_set[[1]] +  -0.002894367 *xyl_set[[1]]^2 + 0.80902937 *yVar$TD+ -0.001691784 *yVar$TD^2 + -0.008893842 *xyl_set[[2]]*xyl_set[[1]]
      dd5$DD5_12= dd5$DD5_12 + 65.34131917 + -3.37385005 *xyl_set[[2]] + 0.014957337 *xyl_set[[2]]^2 + 0.029918066 *xyl_set[[1]] +  -0.004781123 *xyl_set[[1]]^2 + 1.499156184 *yVar$TD+ -0.023509362 *yVar$TD^2 + -0.018515145 *xyl_set[[2]]*xyl_set[[1]]
      dd5$dd5_py12= dd5$dd5_py12 + 65.34131917 + -3.37385005 *xyl_set[[2]] + 0.014957337 *xyl_set[[2]]^2 + 0.029918066 *xyl_set[[1]] +  -0.004781123 *xyl_set[[1]]^2 + 1.499156184 *yVar$TD+ -0.023509362 *yVar$TD^2 + -0.018515145 *xyl_set[[2]]*xyl_set[[1]]
      
      dd5[dd5<0] <- 0; head(dd5)
      dd5$DD5_wt <- rowSums(dd5[,c(1,2,13)]);head(dd5)
      dd5$DD5_sp <- rowSums(dd5[,3:5]);head(dd5)
      dd5$DD5_sm <- rowSums(dd5[,6:8]);head(dd5)
      dd5$DD5_at <- rowSums(dd5[,9:11]);head(dd5)
      dd5$DD5 <- rowSums(dd5[,1:12]);head(dd5)
      dd5 <- round(dd5,0);head(dd5)
      DT0_pdd[row0:row2,names(dd5)[c(1:12,14:18)] := dd5[,c(1:12,14:18)]];DT0_pdd
      rm(dd5);gc()
      
      dd61 <- data.table(DD61_01= (51.8489/(1+exp(-(Tave3[[1]]-	6.46)/1.8843))));head(dd61)
      dd61$DD61_02	=	184.1036	/(1+exp(-(Tave3[[2]]-	9.36)/	1.8712))
      dd61$DD61_03	= ifelse(Tave3[[3]]<=	6.2, 69.7167/(1+exp(-(Tave3[[3]]-	8.12	)/	3.0676	)),	10.0000	*Tave3[[3]]-22.22)
      dd61$DD61_04	=ifelse(	Tave3[[4]]<=	8,	86.8402/(1+exp(-(Tave3[[4	]]-5.27	)/	1.9962	)),	27.7992	*Tave3[[4]]-159.34)
      dd61$DD61_05	=ifelse(	Tave3[[5]]<=	9.3,	260.2962/(1+exp(-(Tave3[[5]]-	9.95	)/	2.4953	)),	27.4567	*Tave3[[5]]+	-143.55)
      dd61$DD61_06	=ifelse(	Tave3[[6]]<=	15.4,	391.1684/(1+exp(-(Tave3[[6]]-	12.67)/	3.2585	)),	30.1523	*Tave3[[6]]+	-190.32)
      dd61$DD61_07	=ifelse(	Tave3[[7]]<=	12,	258.5566	/(1+exp(-(Tave3[[7]]-	10.10	)/2.3152	)),	30.1197	*Tave3[[7]]+	-180.51)
      dd61$DD61_08	=ifelse(	Tave3[[8]]<=	12,	259.1189	/(1+exp(-(Tave3[[8]]-	10.18	)/2.2623	)),	30.3691	*Tave3[[8]]+	-184.69)
      dd61$DD61_09	=ifelse(	Tave3[[9]]<=	10,	214.8702	/(1+exp(-(Tave3[[9]]-	9.15	)/2.2722	)),	27.8971	*Tave3[[9]]+	-155.18)
      dd61$DD61_10	=ifelse(	Tave3[[10]]<=	10,	250.6336	/(1+exp(-(Tave3[[10]]-	10.07)/	3.0549	)),	29.3281	*Tave3[[10]]-167.50)
      dd61$DD61_11	=ifelse(	Tave3[[11]]<=	5,	62.1891	/(1+exp(-(Tave3[[11]]-5.24)/	2.6160)),	10.5930	*Tave3[[11]]	-28.29)
      dd61$DD61_12	=ifelse(	Tave3[[12]]<=	2.5,	9.8083	/(1+exp(-(Tave3[[12]]-1.84)/3.0690)),	5.3790	*Tave3[[12]]-7.65)
      dd61 <- round(dd61,0);head(dd61)
      DT0_pdd[row0:row2,names(dd61)[1:12] := dd61[,1:12]];DT0_pdd
      rm(dd61);gc()
      
      
      ##DD<18 ------------------------------------------
      dd_18 <- data.table(DD_18_01 =ifelse(Tave3[[ 1 ]] >=  11 , 342.1497  / (1 + exp(-(Tave3[[ 1 ]]- 12.78 ) / -3.1966 )), -30.7428  * Tave3[[ 1 ]] + 560 ));head(dd_18)
      dd_18$DD_18_02 =ifelse(Tave3[[ 2 ]] >=  11 , 344.9851  / (1 + exp(-(Tave3[[ 2 ]]- 12.03 ) / -3.2953 )), -28.0059  * Tave3[[ 2 ]] + 500 )
      dd_18$DD_18_03 =ifelse(Tave3[[ 3 ]] >=  11 , 325.8230  / (1 + exp(-(Tave3[[ 3 ]]- 13.10 ) / -2.9954 )), -30.9798  * Tave3[[ 3 ]] + 560 )
      dd_18$DD_18_04 =ifelse(Tave3[[ 4 ]] >=  10 , 325.2590  / (1 + exp(-(Tave3[[ 4 ]]- 12.99 ) / -2.9252 )), -29.8772  * Tave3[[ 4 ]] + 540 )
      dd_18$DD_18_05 =ifelse(Tave3[[ 5 ]] >=  10 , 311.8766  / (1 + exp(-(Tave3[[ 5 ]]- 13.39 ) / -2.7698 )), -30.9508  * Tave3[[ 5 ]] + 558 )
      dd_18$DD_18_06 =ifelse(Tave3[[ 6 ]] >=  12.5 , 220.3158  / (1 + exp(-(Tave3[[ 6 ]]- 14.77 ) / -2.3415 )), -29.9192  * Tave3[[ 6 ]] + 540 )
      dd_18$DD_18_07 =ifelse(Tave3[[ 7 ]] >=  13 , 210.8181  / (1 + exp(-(Tave3[[ 7 ]]- 14.83 ) / -2.0157 )), -31.1228  * Tave3[[ 7 ]] + 560 )
      dd_18$DD_18_08 =ifelse(Tave3[[ 8 ]] >=  14 , 184.0869  / (1 + exp(-(Tave3[[ 8 ]]- 15.45 ) / -2.0060 )), -31.0299  * Tave3[[ 8 ]] + 560 )
      dd_18$DD_18_09 =ifelse(Tave3[[ 9 ]] >=  11 , 298.2082  / (1 + exp(-(Tave3[[ 9 ]]- 13.37 ) / -2.8722 )), -29.9141  * Tave3[[ 9 ]] + 540 )
      dd_18$DD_18_10 =ifelse(Tave3[[ 10 ]] >=  11 , 308.4115  / (1 + exp(-(Tave3[[ 10 ]]- 13.43 ) / -2.8292 )), -31.1362  * Tave3[[ 10 ]] + 560 )
      dd_18$DD_18_11 =ifelse(Tave3[[ 11 ]] >=  11 , 320.2009  / (1 + exp(-(Tave3[[ 11 ]]- 13.01 ) / -3.0693 )), -29.9279  * Tave3[[ 11 ]] + 540 )
      dd_18$DD_18_12 =ifelse(Tave3[[ 12 ]] >=  11 , 353.3562  / (1 + exp(-(Tave3[[ 12 ]]- 12.61 ) / -3.2811 )), -30.9299  * Tave3[[ 12 ]] + 555 )
      dd_18$dd_18_py12 =ifelse(Tave12[[1]] >=  11 , 353.3562  / (1 + exp(-(Tave12[[1]]- 12.61 ) / -3.2811 )), -30.9299  * Tave12[[1]] + 555 );head(dd_18)
      
      dd_18[dd_18<0] <- 0;
      dd_18$DD_18_wt <- rowSums(dd_18[,c(1,2,13)]);
      dd_18$DD_18_sp <- rowSums(dd_18[,3:5]);
      dd_18$DD_18_sm <- rowSums(dd_18[,6:8]);
      dd_18$DD_18_at <- rowSums(dd_18[,9:11]);
      dd_18$DD_18 <- rowSums(dd_18[,1:12]);head(dd_18)
      dd_18 <- round(dd_18,0);head(dd_18)
      DT0_ddnp[row0:row2,names(dd_18)[c(1:12,14:18)] := dd_18[,c(1:12,14:18)]];DT0_ddnp
      rm(dd_18);gc()
      
      ##DD>18 ------------------------------------------
      dd18 <- data.table(DD18_01 =ifelse(Tave3[[ 1]] <=  35 , 252.4909  / (1 + exp(-(Tave3[[ 1]] - 21.48 ) / 2.6795 )), 0.0000  * Tave3[[ 1 ]] + 0 ));head(dd18)
      dd18$DD18_02 =ifelse(Tave3[[ 2]] <=  35 , 154.5938  / (1 + exp(-(Tave3[[ 2]] - 19.81 ) / 2.2733 )), 0.0000  * Tave3[[ 2 ]] + -220 )
      dd18$DD18_03 =ifelse(Tave3[[ 3]] <=  22 , 218.9100  / (1 + exp(-(Tave3[[ 3]] - 20.97 ) / 2.6333 )), 26.5492  * Tave3[[ 3 ]] + -450 )
      dd18$DD18_04 =ifelse(Tave3[[ 4]] <=  23 , 262.9456  / (1 + exp(-(Tave3[[ 4]] - 22.19 ) / 2.9075 )), 28.4619  * Tave3[[ 4 ]] + -500 )
      dd18$DD18_05 =ifelse(Tave3[[ 5]] <=  23 , 270.0578  / (1 + exp(-(Tave3[[ 5]] - 22.16 ) / 2.8396 )), 28.6688  * Tave3[[ 5 ]] + -500 )
      dd18$DD18_06 =ifelse(Tave3[[ 6]] <=  21 , 154.5632  / (1 + exp(-(Tave3[[ 6]] - 19.73 ) / 2.0088 )), 28.8946  * Tave3[[ 6 ]] + -510 )
      dd18$DD18_07 =ifelse(Tave3[[ 7]] <=  22 , 181.7177  / (1 + exp(-(Tave3[[ 7]] - 20.52 ) / 1.8257 )), 30.7090  * Tave3[[ 7 ]] + -550 )
      dd18$DD18_08 =ifelse(Tave3[[ 8]] <=  22 , 190.7156  / (1 + exp(-(Tave3[[ 8]] - 20.63 ) / 1.9841 )), 30.7189  * Tave3[[ 8 ]] + -550 )
      dd18$DD18_09 =ifelse(Tave3[[ 9]] <=  24 , 255.1446  / (1 + exp(-(Tave3[[ 9]] - 21.77 ) / 2.6535 )), 29.2723  * Tave3[[ 9 ]] + -520 )
      dd18$DD18_10 =ifelse(Tave3[[ 10]] <=  23 , 236.9956  / (1 + exp(-(Tave3[[ 10]] - 21.31 ) / 2.5968 )), 28.7088  * Tave3[[ 10 ]] + -500 )
      dd18$DD18_11 =ifelse(Tave3[[ 11]] <=  21 , 159.6297  / (1 + exp(-(Tave3[[ 11]] - 19.52 ) / 2.3321 )), 24.0432  * Tave3[[ 11 ]] + -400 )
      dd18$DD18_12 =ifelse(Tave3[[ 12]] <=  21 , 144.0946  / (1 + exp(-(Tave3[[ 12]] - 18.93 ) / 2.3731 )), 24.2537  * Tave3[[ 12 ]] + -400 )
      dd18$DD18_py12=ifelse(Tave12[[1]] <=  21 , 144.0946  / (1 + exp(-(Tave12[[1]] - 18.93 ) / 2.3731 )), 24.2537  * Tave12[[1]] + -400 );head(dd18)
      
      #spatial adjustments ---
      dd18$DD18_01 =dd18$DD18_01 + 145.6704924 + -13.43358507 *xyl_set[[2]]+ 0.405357077 *xyl_set[[2]]^2+ -0.005227097 *xyl_set[[2]]^3+ 2.46567E-05 *xyl_set[[2]]^4+ -0.110727329 *xyl_set[[1]]+ -0.000659835 *xyl_set[[1]]^2+ 1.662534477 *yVar$TD+ -0.034651808 *yVar$TD^2+ 0.000436466 *yVar$TD^3+ -0.000417514 *xyl_set[[1]]*xyl_set[[2]]+ 0.00448734 *xyl_set[[1]]*yVar$TD+ -0.016121707 *xyl_set[[2]]*yVar$TD+ -8.59015E-05 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD
      dd18$DD18_02 =dd18$DD18_02 + 135.0611504 + -12.16001857 *xyl_set[[2]]+ 0.359336062 *xyl_set[[2]]^2+ -0.004577538 *xyl_set[[2]]^3+ 2.13888E-05 *xyl_set[[2]]^4+ -0.088883611 *xyl_set[[1]]+ -0.000892245 *xyl_set[[1]]^2+ 1.555011969 *yVar$TD+ -0.035390849 *yVar$TD^2+ 0.000457008 *yVar$TD^3+ -0.001711326 *xyl_set[[1]]*xyl_set[[2]]+ 0.003220872 *xyl_set[[1]]*yVar$TD+ -0.014406816 *xyl_set[[2]]*yVar$TD+ -6.0417E-05 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD
      dd18$DD18_03 =dd18$DD18_03 + 86.02322008 + -9.80920465 *xyl_set[[2]]+ 0.286028214 *xyl_set[[2]]^2+ -0.00336782 *xyl_set[[2]]^3+ 1.4981E-05 *xyl_set[[2]]^4+ -0.372993379 *xyl_set[[1]]+ -0.00191843 *xyl_set[[1]]^2+ 3.192734023 *yVar$TD+ -0.041146988 *yVar$TD^2+ 0.000551092 *yVar$TD^3+ 0.00084701 *xyl_set[[1]]*xyl_set[[2]]+ 0.009591184 *xyl_set[[1]]*yVar$TD+ -0.048687383 *xyl_set[[2]]*yVar$TD+ -0.000209062 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD
      dd18$DD18_04 =dd18$DD18_04 + 25.26650509 + -5.140292117 *xyl_set[[2]]+ 0.158088313 *xyl_set[[2]]^2+ -0.001750112 *xyl_set[[2]]^3+ 7.9628E-06 *xyl_set[[2]]^4+ -0.400231953 *xyl_set[[1]]+ -0.001783921 *xyl_set[[1]]^2+ 3.140171201 *yVar$TD+ 0.002769639 *yVar$TD^2+ 3.44861E-05 *yVar$TD^3+ 0.002523367 *xyl_set[[1]]*xyl_set[[2]]+ 0.012335185 *xyl_set[[1]]*yVar$TD+ -0.071451469 *xyl_set[[2]]*yVar$TD+ -0.000282954 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD
      dd18$DD18_05 =dd18$DD18_05 + 211.8587556 + -21.0942444 *xyl_set[[2]]+ 0.662425055 *xyl_set[[2]]^2+ -0.008677058 *xyl_set[[2]]^3+ 4.2388E-05 *xyl_set[[2]]^4+ -0.196878146 *xyl_set[[1]]+ -0.000428519 *xyl_set[[1]]^2+ 3.27379518 *yVar$TD+ 0.013835765 *yVar$TD^2+ -5.1475E-06 *yVar$TD^3+ 0.001851708 *xyl_set[[1]]*xyl_set[[2]]+ 0.013802948 *xyl_set[[1]]*yVar$TD+ -0.071078166 *xyl_set[[2]]*yVar$TD+ -0.000250307 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD
      dd18$DD18_10 =dd18$DD18_10 + 65.57254896 + -6.047086844 *xyl_set[[2]]+ 0.140233853 *xyl_set[[2]]^2+ -0.001231275 *xyl_set[[2]]^3+ 4.4462E-06 *xyl_set[[2]]^4+ -0.055221198 *xyl_set[[1]]+ -0.002208941 *xyl_set[[1]]^2+ 2.293734141 *yVar$TD+ -0.013108002 *yVar$TD^2+ 0.000238296 *yVar$TD^3+ -0.005876538 *xyl_set[[1]]*xyl_set[[2]]+ -0.004059974 *xyl_set[[1]]*yVar$TD+ -0.048245454 *xyl_set[[2]]*yVar$TD+ 2.6668E-05 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD
      dd18$DD18_11 =dd18$DD18_11 + 191.1500248 + -20.0095113 *xyl_set[[2]]+ 0.598904826 *xyl_set[[2]]^2+ -0.007458949 *xyl_set[[2]]^3+ 3.43629E-05 *xyl_set[[2]]^4+ -0.493033407 *xyl_set[[1]]+ -0.001991277 *xyl_set[[1]]^2+ 4.664022754 *yVar$TD+ -0.079247685 *yVar$TD^2+ 0.001014898 *yVar$TD^3+ 0.002024583 *xyl_set[[1]]*xyl_set[[2]]+ 0.015266276 *xyl_set[[1]]*yVar$TD+ -0.055272548 *xyl_set[[2]]*yVar$TD+ -0.000290957 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD
      dd18$DD18_12 =dd18$DD18_12 + 129.6810431 + -12.39710088 *xyl_set[[2]]+ 0.370770764 *xyl_set[[2]]^2+ -0.004744748 *xyl_set[[2]]^3+ 2.22776E-05 *xyl_set[[2]]^4+ -0.134456142 *xyl_set[[1]]+ -0.000887304 *xyl_set[[1]]^2+ 2.131816236 *yVar$TD+ -0.049138542 *yVar$TD^2+ 0.000627762 *yVar$TD^3+ -0.001036013 *xyl_set[[1]]*xyl_set[[2]]+ 0.004878242 *xyl_set[[1]]*yVar$TD+ -0.018075939 *xyl_set[[2]]*yVar$TD+ -8.47361E-05 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD
      dd18$DD18_py12 =dd18$DD18_py12 + 129.6810431 + -12.39710088 *xyl_set[[2]]+ 0.370770764 *xyl_set[[2]]^2+ -0.004744748 *xyl_set[[2]]^3+ 2.22776E-05 *xyl_set[[2]]^4+ -0.134456142 *xyl_set[[1]]+ -0.000887304 *xyl_set[[1]]^2+ 2.131816236 *yVar$TD+ -0.049138542 *yVar$TD^2+ 0.000627762 *yVar$TD^3+ -0.001036013 *xyl_set[[1]]*xyl_set[[2]]+ 0.004878242 *xyl_set[[1]]*yVar$TD+ -0.018075939 *xyl_set[[2]]*yVar$TD+ -8.47361E-05 *xyl_set[[2]]*xyl_set[[1]]*yVar$TD;head(dd18)
      
      dd18[dd18<1] <- 0;
      dd18$DD18_wt <- rowSums(dd18[,c(1,2,13)]);
      dd18$DD18_sp <- rowSums(dd18[,3:5]);
      dd18$DD18_sm <- rowSums(dd18[,6:8]);
      dd18$DD18_at <- rowSums(dd18[,9:11]);
      dd18$DD18 <- rowSums(dd18[,1:12]);
      dd18 <- round(dd18,0);head(dd18)
      DT0_ddnp[row0:row2,names(dd18)[c(1:12,14:18)] := dd18[,c(1:12,14:18)]];DT0_ddnp
      rm(dd18);gc()
      
      #NFFD ------------------------------------------------------------------
      nffd_wna <- data.table(NFFD01  = 34.742  / (1 + exp(-(Tmin3[[ 1 ]] - 0.86 ) / 3.1104 )));head(nffd_wna)
      nffd_wna$NFFD02 = 30.3854 / (1 + exp(-(Tmin3[[2]] - 0.63) / 2.7761))
      nffd_wna$NFFD03 = 32.5385 / (1 + exp(-(Tmin3[[3]] - 0.45) / 2.6121))
      nffd_wna$NFFD04 = 30.6635 / (1 + exp(-(Tmin3[[4]] - 0.36) / 2.3966))
      nffd_wna$NFFD05 = 31.2419 / (1 + exp(-(Tmin3[[5]] - 0.47) / 2.2458))
      nffd_wna$NFFD06 = 30.0486 / (1 + exp(-(Tmin3[[6]] - 0.4) / 1.9178))
      nffd_wna$NFFD07 = 30.9662 / (1 + exp(-(Tmin3[[7]] - 0.67) / 1.4936))
      nffd_wna$NFFD08 = 30.9824 / (1 + exp(-(Tmin3[[8]] - 0.55) / 1.7314))
      nffd_wna$NFFD09 = 30.1372 / (1 + exp(-(Tmin3[[9]] - 0.33) / 2.46))
      nffd_wna$NFFD10 = 31.5636 / (1 + exp(-(Tmin3[[10]] - 0.34) / 2.548))
      nffd_wna$NFFD11 = 31.6777 / (1 + exp(-(Tmin3[[11]] - 0.45) / 2.8848))
      nffd_wna$NFFD12 = 34.409 / (1 + exp(-(Tmin3[[12]] - 0.78) / 2.9289))
      nffd_wna$NFFD_py12 = 34.409 / (1 + exp(-(Tmin12[[1]] - 0.78) / 2.9289));head(nffd_wna)
      
      nffd_NA <- data.table(NFFD01  = 31.9203  / (1 + exp(-(Tmin3[[ 1 ]] - 0.96 ) / 3.8184 )));head(nffd_NA)
      nffd_NA$NFFD02  = 29.4221  / (1 + exp(-(Tmin3[[ 2 ]] - 3.80 ) / 0.6279 ))
      nffd_NA$NFFD03  = 31.9966  / (1 + exp(-(Tmin3[[ 3 ]] - 3.60 ) / 0.7681 ))
      nffd_NA$NFFD04  = 30.4145  / (1 + exp(-(Tmin3[[ 4 ]] - 3.13 ) / 0.8780 ))
      nffd_NA$NFFD05  = 31.2379  / (1 + exp(-(Tmin3[[ 5 ]] - 2.78 ) / 0.7785 ))
      nffd_NA$NFFD06  = 30.0053  / (1 + exp(-(Tmin3[[ 6 ]] - 2.24 ) / 0.3261 ))
      nffd_NA$NFFD07  = 30.9517  / (1 + exp(-(Tmin3[[ 7 ]] - 1.23 ) / 0.1754 ))
      nffd_NA$NFFD08  = 30.9461  / (1 + exp(-(Tmin3[[ 8 ]] - 1.71 ) / 0.2018 ))
      nffd_NA$NFFD09  = 30.1120  / (1 + exp(-(Tmin3[[ 9 ]] - 2.72 ) / 0.4759 ))
      nffd_NA$NFFD10  = 31.5968  / (1 + exp(-(Tmin3[[ 10 ]] - 3.23 ) / 0.8385 ))
      nffd_NA$NFFD11  = 30.5354  / (1 + exp(-(Tmin3[[ 11 ]] - 3.47 ) / 0.7373 ))
      nffd_NA$NFFD12  = 31.3262  / (1 + exp(-(Tmin3[[ 12 ]] - 3.63 ) / 0.7792 ))
      nffd_NA$NFFD_py12 = 31.3262  / (1 + exp(-(Tmin12 - 3.63 ) / 0.7792 ));head(nffd_NA)
      
      #Combine wna and na
      Wt_wna = abs(xyl_set$x - (-95)) / 10; head(Wt_wna)
      Wt_na = 1 - Wt_wna; head(Wt_na)
      nffd <- data.table(NFFD01 = ifelse(xyl_set$x <= -105, nffd_wna$NFFD01, nffd_NA$NFFD01));head(nffd)
      nffd$NFFD02 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD02, nffd_NA$NFFD02);head(nffd)
      nffd$NFFD03 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD03, nffd_NA$NFFD03);head(nffd)
      nffd$NFFD04 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD04, nffd_NA$NFFD04);head(nffd)
      nffd$NFFD05 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD05, nffd_NA$NFFD05);head(nffd)
      nffd$NFFD06 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD06, nffd_NA$NFFD06);head(nffd)
      nffd$NFFD07 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD07, nffd_NA$NFFD07);head(nffd)
      nffd$NFFD08 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD08, nffd_NA$NFFD08);head(nffd)
      nffd$NFFD09 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD09, nffd_NA$NFFD09);head(nffd)
      nffd$NFFD10 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD10, nffd_NA$NFFD10);head(nffd)
      nffd$NFFD11 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD11, nffd_NA$NFFD11);head(nffd)
      nffd$NFFD12 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD12, nffd_NA$NFFD12);head(nffd)
      nffd$NFFD_py12 <- ifelse(xyl_set$x <= -105, nffd_wna$NFFD_py12, nffd_NA$NFFD_py12);head(nffd)
      
      nffd$NFFD01 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD01*Wt_wna+nffd_NA$NFFD01 * Wt_na,nffd$NFFD01);head(nffd)
      nffd$NFFD02 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD02*Wt_wna+nffd_NA$NFFD02 * Wt_na,nffd$NFFD02);head(nffd)
      nffd$NFFD03 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD03*Wt_wna+nffd_NA$NFFD03 * Wt_na,nffd$NFFD03);head(nffd)
      nffd$NFFD04 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD04*Wt_wna+nffd_NA$NFFD04 * Wt_na,nffd$NFFD04);head(nffd)
      nffd$NFFD05 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD05*Wt_wna+nffd_NA$NFFD05 * Wt_na,nffd$NFFD05);head(nffd)
      nffd$NFFD06 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD06*Wt_wna+nffd_NA$NFFD06 * Wt_na,nffd$NFFD06);head(nffd)
      nffd$NFFD07 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD07*Wt_wna+nffd_NA$NFFD07 * Wt_na,nffd$NFFD07);head(nffd)
      nffd$NFFD08 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD08*Wt_wna+nffd_NA$NFFD08 * Wt_na,nffd$NFFD08);head(nffd)
      nffd$NFFD09 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD09*Wt_wna+nffd_NA$NFFD09 * Wt_na,nffd$NFFD09);head(nffd)
      nffd$NFFD10 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD10*Wt_wna+nffd_NA$NFFD10 * Wt_na,nffd$NFFD10);head(nffd)
      nffd$NFFD11 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD11*Wt_wna+nffd_NA$NFFD11 * Wt_na,nffd$NFFD11);head(nffd)
      nffd$NFFD12 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD12*Wt_wna+nffd_NA$NFFD12 * Wt_na,nffd$NFFD12);head(nffd)
      nffd$NFFD_py12 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, nffd_wna$NFFD_py12*Wt_wna+nffd_NA$NFFD_py12 * Wt_na,nffd$NFFD_py12);head(nffd)
      
      nffd[nffd<0] <- 0
      nffd$NFFD01 <- ifelse(nffd$NFFD01>31,31,nffd$NFFD01);head(nffd)
      nffd$NFFD02 <- ifelse(nffd$NFFD02>28,28,nffd$NFFD02);head(nffd)
      nffd$NFFD03 <- ifelse(nffd$NFFD03>31,31,nffd$NFFD03);head(nffd)
      nffd$NFFD04 <- ifelse(nffd$NFFD04>30,30,nffd$NFFD04);head(nffd)
      nffd$NFFD05 <- ifelse(nffd$NFFD05>31,31,nffd$NFFD05);head(nffd)
      nffd$NFFD06 <- ifelse(nffd$NFFD06>30,30,nffd$NFFD06);head(nffd)
      nffd$NFFD07 <- ifelse(nffd$NFFD07>31,31,nffd$NFFD07);head(nffd)
      nffd$NFFD08 <- ifelse(nffd$NFFD08>31,31,nffd$NFFD08);head(nffd)
      nffd$NFFD09 <- ifelse(nffd$NFFD09>30,30,nffd$NFFD09);head(nffd)
      nffd$NFFD10 <- ifelse(nffd$NFFD10>31,31,nffd$NFFD10);head(nffd)
      nffd$NFFD11 <- ifelse(nffd$NFFD11>30,30,nffd$NFFD11);head(nffd)
      nffd$NFFD12 <- ifelse(nffd$NFFD12>31,31,nffd$NFFD12);head(nffd)
      nffd$NFFD_py12 <- ifelse(nffd$NFFD_py12>31,31,nffd$NFFD_py12);head(nffd)
      
      nffd$NFFD_wt <- rowSums(nffd[,c(1,2,13)]);nffd
      nffd$NFFD_sp <- rowSums(nffd[,3:5]);nffd
      nffd$NFFD_sm <- rowSums(nffd[,6:8]);nffd
      nffd$NFFD_at <- rowSums(nffd[,9:11]);nffd
      nffd$NFFD <- rowSums(nffd[,1:12]);nffd
      nffd <- round(nffd,0);head(nffd)
      DT0_ddnp[row0:row2,names(nffd)[c(1:12,14:18)] := nffd[,c(1:12,14:18)]];DT0_ddnp
      rm(nffd_NA,nffd_wna); gc()
      
      #bFFP--------------->>
      txn04 = Tmax3[[4]] - Tmin3[[4]];head(txn04)
      bFFP_wna = 124.9495 + Tmin3[[4]] * -1.7581+Tmin3[[5]] * -11.87934+Tmin3[[6]]*2.09433+Tmin3[[4]]^2 * -0.3746 + Tmin3[[4]]^3 * 0.01482 + Tmin3[[5]] ^ 3 * 0.06751 + Tmin3[[4]]^4 * 0.00123+Tmin3[[5]]^4* -0.00266+txn04*5.21934+txn04^2* -0.16101 + nffd$NFFD^3* -0.00000719 + nffd$NFFD^4*0.00000005976953+nffd$NFFD^5* -0.00000000012266;head(bFFP_wna)
      bFFP_na  = 352.1358994+ -0.021715653*Tmin3[[4]]^2+ -3.542187618*Tmin3[[6]]+0.020359471*Tmin3[[6]]^2 + -4.897998097*yVar$TD+0.033521327*yVar$TD^2+ -2.164862277*nffd$NFFD+0.006767633*nffd$NFFD^2+ -0.00000929*nffd$NFFD^3+0.043516586*(yVar$TD*nffd$NFFD)+ -0.00000253*(yVar$TD*nffd$NFFD)^2;head(bFFP_na)
      
      #Combine wna and na
      ffp <- data.table(bFFP=ifelse(xyl_set$x <= -105, bFFP_wna,bFFP_na));head(ffp)
      ffp$bFFP <- ifelse((xyl_set$x > -105 & xyl_set$x < -95), bFFP_wna * Wt_wna + bFFP_na * Wt_na,ffp$bFFP);head(ffp) #WT are calculated in NFFD
      ffp[bFFP< 0] <- 0; ffp[bFFP > 365] <- 365; head(ffp)
      
      #eFFP -------------->>
      eFFP_wna=231.65777+Tmin3[[9]]*8.87656+Tmin3[[6]]^2* -0.05996+Tmin3[[7]]^2* -0.0751+Tmin3[[10]]^2*0.20123+Tmin3[[9]]^3* -0.026+Tmin3[[7]]^4*0.00009435+Tmin3[[9]]^4*0.00067816+Tmin3[[10]]^4* -0.00029319+nffd$NFFD^5*0.00000000000794
      eFFP_na <- 243.7752209 + 4.134210825 * Tmin3[[9]] + -0.162876448 * Tmin3[[9]]^2 + 1.248649021 * Tmin3[[10]] + 0.145073612 * Tmin3[[10]]^2 + 0.004319892 * Tmin3[[11]] + -0.005753127 * Tmin3[[11]]^2 + -0.062964718 * nffd$NFFD + 0.000399177 * nffd$NFFD^2;
      ffp$eFFP <- ifelse(xyl_set$x <= -105, eFFP_wna,eFFP_na);head(ffp);length(ffp)
      ffp$eFFP <- ifelse(xyl_set$x >  -105 & xyl_set$x < -95, eFFP_wna * Wt_wna + eFFP_na * Wt_na,ffp$eFFP);head(ffp) #WT are calculated in NFFD
      
      #FFP -----------
      #ffp <- data.table(bFFP,eFFP);head(ffp)
      ffp$FFP <- ffp$eFFP - ffp$bFFP;head(ffp)
      ffp[ffp<0] <- 0; head(ffp)
      ffp[ffp>365] <- 365;head(ffp)
      ffp <- round(ffp,0);head(ffp)
      DT0_ecr[row0:row2,names(ffp) := ffp];DT0_ecr
      rm(txn04,bFFP_na,bFFP_wna,eFFP_wna,eFFP_na,ffp);gc()
      
      ##EMT and EXT --------------------
      #EMT
      #EMT for WNA
      EMT_wna = -18.09995 + Tmin3[[1]] * 2.14095 + Tmin3[[1]]^2 * 0.06836 + Tmin3[[12]]^2* -0.04771 + yVar$TD^2 * 0.00306
      if(prd_class=='His_Ann'){
        TD_nrm = (Tmax3_cruNrm[[7]] + Tmin3_cruNrm[[7]]) / 2 - (Tmax3_cruNrm[[1]] + Tmin3_cruNrm[[1]]) / 2
        EMT_wna = -18.09995 + Tmin3_cruNrm[[1]] * 2.14095 + Tmin3_cruNrm[[1]]^2 * 0.06836 + Tmin3_cruNrm[[12]]^2 * -0.04771 + TD_nrm ^ 2 * 0.00306
      }
      
      #EMT for the rest of NA
      Tminx = apply(Tmin3, 1, FUN = min);head(Tminx)
      EMT_na = -23.02164 + 0.77908 * Tmin3[[1]] + 0.67048 * Tmin3[[12]] + 0.01075 * Tminx^2 + 0.11565 * yVar$TD
      if(prd_class =='His_Ann'){
        TD_nrm = (Tmax3_cruNrm[[7]] + Tmin3_cruNrm[[7]]) / 2 - (Tmax3_cruNrm[[1]] + Tmin3_cruNrm[[1]]) / 2
        Tminx_nrm = apply(Tmin3_cruNrm, 1, FUN = min);head(Tminx_nrm)
        EMT_na = -23.02164 + 0.77908 * Tmin3_cruNrm[[1]]+0.67048*Tmin3_cruNrm[[12]] + 0.01075 * Tminx_nrm^2 + 0.11565 * TD_nrm
      }
      
      #combine
      EMXT = data.table(EMT=ifelse(xyl_set$x <= -105, EMT_wna,EMT_na));head(EMXT)
      EMXT$EMT <- ifelse(xyl_set$x > -105 & xyl_set$x < -95, EMT_wna * Wt_wna + EMT_na * Wt_na,EMXT[[1]]);head(EMXT) #WT are calculated in NFFD
      rm(EMT_na,EMT_wna);gc()
      
      #EXT --------------
      #EXT_wna
      EXT_wna  = 16.785 + 0.5252 * Tmax3[[7]] - 0.0294 * Tmax3[[7]]^2 + 0.4391 * Tmax3[[8]] + 0.028 * Tmax3[[8]]^2 - 0.3191 * yVar$TD + 0.00813 * yVar$TD^2
      if(prd_class =='His_Ann'){
        EXT_wna = 16.785 + 0.5252 * Tmax3_cruNrm[[7]] - 0.0294 * Tmax3_cruNrm[[7]] ^ 2 + 0.4391 * Tmax3_cruNrm[[8]] + 0.028 * Tmax3_cruNrm[[8]] ^ 2 - 0.3191 * TD_nrm + 0.00813 * TD_nrm ^ 2
      }
      
      #EXT_na
      Tmaxx = apply(Tmax3, 1, FUN = max);head(Tmaxx)
      EXT_na = 10.64245 + -1.92005 * Tmax3[[7]] + 0.04816 * Tmax3[[7]]^2 + 2.51176 * Tmax3[[8]] + -0.03088 * Tmax3[[8]]^2 + -0.01311 * Tmaxx^2 + 0.33167 * yVar$TD + -0.001 * yVar$TD^2;
      if(prd_class =='His_Ann'){
        Tmaxx_nrm = apply(Tmax3_cruNrm, 1, FUN = max);head(Tmaxx_nrm)
        EXT_na = 10.64245 + -1.92005 * Tmax3_cruNrm[[7]] + 0.04816 * Tmax3_cruNrm[[7]]^2 + 2.51176 * Tmax3_cruNrm[[8]]+ -0.03088 * Tmax3_cruNrm[[8]]^2 + -0.01311 * Tmaxx_nrm ^ 2 + 0.33167 * TD_nrm + -0.001 * TD_nrm ^ 2
      }
      
      #combine
      EMXT$EXT <- ifelse(xyl_set$x <= -105, EXT_wna,EXT_na);head(EMXT)
      EMXT$EXT <- ifelse(xyl_set$x > -105 & xyl_set$x < -95, EXT_wna * Wt_wna + EXT_na * Wt_na,EMXT$EXT);head(EMXT) #WT are calculated in NFFD
      DT1_ttr[row0:row2, names(EMXT) := round(EMXT*10,0)];DT1_ttr
      rm(EXT_na,EXT_wna,EMXT);gc();ls()
      
      #RH ---------------------------------------------------------
      RH <- calculate_RH(Tmax3, Tmin3);head(RH)
      RH12py <- calculate_RH(Tmax12, Tmin12);head(RH12py)
      RH$RH_wt = (RH[,1]+RH[,2]+RH12py)/3;head(RH)
      RH$RH_sp = rowMeans(RH[,c(3:5)]);head(RH)
      RH$RH_sm = rowMeans(RH[,c(6:8)]);head(RH)
      RH$RH_at = rowMeans(RH[,c(9:11)]);head(RH)
      RH$RH <- (rowSums(RH[,c(1:11)]) + RH12py)/12;head(RH)
      RH[RH>100] <- 100
      RH <- round(RH,0);head(RH)
      rh <- data.table(RH)
      DT0_ecr[row0:row2, names(rh) := rh];DT0_ecr
      rm(RH,rh);gc()
      
      
      ##PAS ----------------------------------------------------------------------
      PAS_wna <- data.table(PAS01 = 1 / (1 + exp(-(Tave3[[1]] - -2.9901) / -2.50353)) * ppt3[[1]])
      PAS_wna$PAS02 = 1 / (1 + exp(-(Tave3[[2]] - -1.3948) / -2.0004)) * ppt3[[2]]
      PAS_wna$PAS03 = 1 / (1 + exp(-(Tave3[[3]] - 0.5473) / -1.5719)) * ppt3[[3]]
      PAS_wna$PAS04 = 1 / (1 + exp(-(Tave3[[4]] - 2.0928) / -1.6527)) * ppt3[[4]]
      PAS_wna$PAS05 = 0.8 / (1 + exp(-(Tave3[[5]] - 4.078) / -1.7428)) * ppt3[[5]]
      PAS_wna$PAS06 = 0.8 / (1 + exp(-(Tave3[[6]] - 4.078) / -1.7428)) * ppt3[[6]]
      PAS_wna$PAS07 = 0.9 / (1 + exp(-(Tave3[[7]] - 2.8) / -2.3948)) * ppt3[[7]]
      PAS_wna$PAS08 = 1 / (1 + exp(-(Tave3[[8]] - 1.4927) / -2.8948)) * ppt3[[8]]
      PAS_wna$PAS09 = 1 / (1 + exp(-(Tave3[[9]] - 1.4927) / -2.8948)) * ppt3[[9]]
      PAS_wna$PAS10 = 1 / (1 + exp(-(Tave3[[10]] - 0.8099) / -1.6612)) * ppt3[[10]]
      PAS_wna$PAS11 = 1 / (1 + exp(-(Tave3[[11]] - -1.5627) / -2.4907)) * ppt3[[11]]
      PAS_wna$PAS12 = 1 / (1 + exp(-(Tave3[[12]] - -2.5909) / -2.2108)) * ppt3[[12]]; head(PAS_wna)
      
      PAS_na <- data.table(PAS01 =  1/ (1 + exp(-(Tave3[[ 1 ]] - -4.16 ) / -2.5114 ))  *ppt3[[ 1 ]]);head(PAS_na)
      PAS_na$PAS02 =  1/ (1 + exp(-(Tave3[[ 2 ]] - -2.70 ) / -1.7031 ))  *ppt3[[ 2 ]]
      PAS_na$PAS03 =  1/ (1 + exp(-(Tave3[[ 3 ]] - -1.79 ) / -1.2583 ))  *ppt3[[ 3 ]]
      PAS_na$PAS04 =  1/ (1 + exp(-(Tave3[[ 4 ]] - 1.77 ) / -1.4152 ))  *ppt3[[ 4 ]]
      PAS_na$PAS05 =  1/ (1 + exp(-(Tave3[[ 5 ]] - 1.44 ) / -2.2797 ))  *ppt3[[ 5 ]]
      PAS_na$PAS06 =  1/ (1 + exp(-(Tave3[[ 6 ]] - 1.44 ) / -2.2797 ))  *ppt3[[ 6 ]]
      PAS_na$PAS07 =  1/ (1 + exp(-(Tave3[[ 7 ]] - 2.3  ) / -2.1    ))  *ppt3[[ 7 ]]
      PAS_na$PAS08 =  1/ (1 + exp(-(Tave3[[ 8 ]] - 3.2  ) / -1.9808 ))  *ppt3[[ 8 ]]
      PAS_na$PAS09 =  1/ (1 + exp(-(Tave3[[ 9 ]] - 3.20 ) / -1.9808 ))  *ppt3[[ 9 ]]
      PAS_na$PAS10 =  1/ (1 + exp(-(Tave3[[ 10 ]] - 2.35 ) / -1.4464 ))  *ppt3[[ 10 ]]
      PAS_na$PAS11 =  1/ (1 + exp(-(Tave3[[ 11 ]] - -1.67 ) / -1.4617 ))  *ppt3[[ 11 ]]
      PAS_na$PAS12 =  1/ (1 + exp(-(Tave3[[ 12 ]] - -3.01 ) / -1.5327 ))  *ppt3[[ 12 ]]; head(PAS_na)
      
      #Combine wna and na
      pas <- data.table(PAS01 = ifelse(xyl_set$x <= -105, PAS_wna$PAS01, PAS_na$PAS01));head(pas)
      pas$PAS02 = ifelse(xyl_set$x <= -105, PAS_wna$PAS02, PAS_na$PAS02);head(pas)
      pas$PAS03 = ifelse(xyl_set$x <= -105, PAS_wna$PAS03, PAS_na$PAS03);head(pas)
      pas$PAS04 = ifelse(xyl_set$x <= -105, PAS_wna$PAS04, PAS_na$PAS04);head(pas)
      pas$PAS05 = ifelse(xyl_set$x <= -105, PAS_wna$PAS05, PAS_na$PAS05);head(pas)
      pas$PAS06 = ifelse(xyl_set$x <= -105, PAS_wna$PAS06, PAS_na$PAS06);head(pas)
      pas$PAS07 = ifelse(xyl_set$x <= -105, PAS_wna$PAS07, PAS_na$PAS07);head(pas)
      pas$PAS08 = ifelse(xyl_set$x <= -105, PAS_wna$PAS08, PAS_na$PAS08);head(pas)
      pas$PAS09 = ifelse(xyl_set$x <= -105, PAS_wna$PAS09, PAS_na$PAS09);head(pas)
      pas$PAS10 = ifelse(xyl_set$x <= -105, PAS_wna$PAS10, PAS_na$PAS10);head(pas)
      pas$PAS11 = ifelse(xyl_set$x <= -105, PAS_wna$PAS11, PAS_na$PAS11);head(pas)
      pas$PAS12 = ifelse(xyl_set$x <= -105, PAS_wna$PAS12, PAS_na$PAS12);head(pas)
      
      pas$PAS01 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS01*Wt_wna+PAS_na$PAS01 * Wt_na,pas$PAS01);head(pas)
      pas$PAS02 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS02*Wt_wna+PAS_na$PAS02 * Wt_na,pas$PAS02);head(pas)
      pas$PAS03 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS03*Wt_wna+PAS_na$PAS03 * Wt_na,pas$PAS03);head(pas)
      pas$PAS04 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS04*Wt_wna+PAS_na$PAS04 * Wt_na,pas$PAS04);head(pas)
      pas$PAS05 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS05*Wt_wna+PAS_na$PAS05 * Wt_na,pas$PAS05);head(pas)
      pas$PAS06 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS06*Wt_wna+PAS_na$PAS06 * Wt_na,pas$PAS06);head(pas)
      pas$PAS07 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS07*Wt_wna+PAS_na$PAS07 * Wt_na,pas$PAS07);head(pas)
      pas$PAS08 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS08*Wt_wna+PAS_na$PAS08 * Wt_na,pas$PAS08);head(pas)
      pas$PAS09 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS09*Wt_wna+PAS_na$PAS09 * Wt_na,pas$PAS09);head(pas)
      pas$PAS10 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS10*Wt_wna+PAS_na$PAS10 * Wt_na,pas$PAS10);head(pas)
      pas$PAS11 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS11*Wt_wna+PAS_na$PAS11 * Wt_na,pas$PAS11);head(pas)
      pas$PAS12 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna$PAS12*Wt_wna+PAS_na$PAS12 * Wt_na,pas$PAS12);head(pas)
      
      pas[pas<0] <- 0; pas;head(pas)
      pas$PAS_wt = rowSums(pas[,c(1,2,12)]);
      pas$PAS_sp <- rowSums(pas[,3:5]);
      pas$PAS_sm <- rowSums(pas[,6:8]);
      pas$PAS_at <- rowSums(pas[,9:11]);
      pas$PAS <- rowSums(pas[,1:12]);head(pas)
      rm(PAS_wna,PAS_na);gc()
      
      #Annual considering previous year PAS
      if(prd_class=='His_Ann' & yr > 1901){
        PAS_wna_py <- data.table(PAS08 = 1 /(1 + exp(-(Tave3_py[[8]] - 1.4927) / -2.8948)) * ppt3_py[[8]])
        PAS_wna_py$PAS09 = 1 / (1 + exp(-(Tave3_py[[9]] - 1.4927) / -2.8948)) * ppt3_py[[9]]
        PAS_wna_py$PAS10 = 1 / (1 + exp(-(Tave3_py[[10]] - 0.8099) / -1.6612)) * ppt3_py[[10]]
        PAS_wna_py$PAS11 = 1 / (1 + exp(-(Tave3_py[[11]] - -1.5627) / -2.4907)) * ppt3_py[[11]]
        PAS_wna_py$PAS12 = 1 / (1 + exp(-(Tave3_py[[12]] - -2.5909) / -2.2108)) * ppt3_py[[12]]
        PAS_na_py <- data.table(PAS08 = 1 / (1 + exp(-(Tave3_py[[8]] - 3.2) / -1.9808)) * ppt3_py[[8]])
        PAS_na_py$PAS09 = 1 / (1 + exp(-(Tave3_py[[9]] - 3.2) / -1.9808)) * ppt3_py[[9]]
        PAS_na_py$PAS10 = 1 / (1 + exp(-(Tave3_py[[10]] - 2.35) / -1.4464)) * ppt3_py[[10]]
        PAS_na_py$PAS11 = 1 / (1 + exp(-(Tave3_py[[11]] - -1.67) / -1.4617)) * ppt3_py[[11]]
        PAS_na_py$PAS12 = 1 / (1 + exp(-(Tave3_py[[12]] - -3.01) / -1.5327)) * ppt3_py[[12]]
        #combine
        pas_py <- data.table(PAS08 = ifelse(xyl_set$x <= -105, PAS_wna_py$PAS08, PAS_na_py$PAS08));head(pas_py)
        pas_py$PAS09 = ifelse(xyl_set$x <= -105, PAS_wna_py$PAS09, PAS_na_py$PAS09);head(pas_py)
        pas_py$PAS10 = ifelse(xyl_set$x <= -105, PAS_wna_py$PAS10, PAS_na_py$PAS10);head(pas_py)
        pas_py$PAS11 = ifelse(xyl_set$x <= -105, PAS_wna_py$PAS11, PAS_na_py$PAS11);head(pas_py)
        pas_py$PAS12 = ifelse(xyl_set$x <= -105, PAS_wna_py$PAS12, PAS_na_py$PAS12);head(pas_py)
        
        pas_py$PAS08 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna_py$PAS08*Wt_wna+PAS_na_py$PAS08 * Wt_na,pas_py$PAS08);head(pas_py)
        pas_py$PAS09 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna_py$PAS09*Wt_wna+PAS_na_py$PAS09 * Wt_na,pas_py$PAS09);head(pas_py)
        pas_py$PAS10 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna_py$PAS10*Wt_wna+PAS_na_py$PAS10 * Wt_na,pas_py$PAS10);head(pas_py)
        pas_py$PAS11 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna_py$PAS11*Wt_wna+PAS_na_py$PAS11 * Wt_na,pas_py$PAS11);head(pas_py)
        pas_py$PAS12 = ifelse(xyl_set$x > -105 & xyl_set$x < -95, PAS_wna_py$PAS12*Wt_wna+PAS_na_py$PAS12 * Wt_na,pas_py$PAS12);head(pas_py)
        
        pas_py[pas_py<0] <- 0; head(pas_py)
        pas$PAS_wt = rowSums(cbind(pas[,c(1,2)],pas_py[,5]));
        pas$PAS_sp <- rowSums(pas[,3:5]);
        pas$PAS_sm <- rowSums(pas[,6:8]);
        pas$PAS_at <- rowSums(pas[,9:11]);
        pas$PAS <- rowSums(cbind(pas_py,pas[,1:7]));head(pas)
        rm(pas_py,PAS_wna_py,PAS_na_py);gc();ls()
      }
      pas = round(pas,0);head(pas)
      DT0_ddnp[row0:row2,names(pas) := pas];DT0_ddnp
      rm(pas);gc()
      
      #Eref ---------------------------------------------
      ecd <- calculate_Eref_CMD(xyl_set, Tmax3, Tmin3, Tave3, ppt3) ;head(ecd)
      Eref <- ecd$Eref;head(Eref)
      ErefR <- ecd$ErefR;head(ErefR)
      CMD <- ecd$CMD;head(CMD)
      
      #previous year --
      ecd12 <- calculate_Eref_CMD(xyl_set, Tmax12, Tmin12, Tave12, PPT12) ;head(ecd12)
      Eref12py <- ecd12$Eref;head(Eref12py)
      ErefR12py <- ecd12$ErefR;head(ErefR12py)
      cmd12py <- ecd12$CMD;head(cmd12py)
      
      #seasonal var -------------
      Eref$Eref_wt <- (Eref[,1]+Eref[,2]+ Eref12py)/3;Eref
      Eref$Eref_sp <- rowSums(Eref[,3:5]);Eref
      Eref$Eref_sm <- rowSums(Eref[,6:8]);Eref
      Eref$Eref_at <- rowSums(Eref[,9:11]);Eref
      Eref$Eref <- rowSums(Eref[,1:12]);Eref
      Eref <- round(Eref,0);Eref
      
      CMD$CMD_wt <- (CMD[,1]+CMD[,2]+cmd12py)/3;CMD
      CMD$CMD_sp <- rowSums(CMD[,3:5]);CMD
      CMD$CMD_sm <- rowSums(CMD[,6:8]);CMD
      CMD$CMD_at <- rowSums(CMD[,9:11]);CMD
      CMD$CMD <- rowSums(CMD[,1:12]);CMD
      CMD <- round(CMD,0);CMD
      
      # add to DT0_ecr -----
      eref2 <- data.table(Eref)
      cmd2 <- data.table(CMD)
      erefr2 <- data.table(ErefR)
      DT0_ecr[row0:row2,names(Eref) := eref2];head(DT0_ecr)
      DT0_ecr[row0:row2,names(cmd2) := cmd2];DT0_ecr;
      DT0_ecr[row0:row2,names(ErefR) := round(erefr2,0)];DT0_ecr;
      rm(Eref,eref2,ErefR,erefr2,CMD,cmd2);gc()
      
      
      #DD1040 --------------------------------------
      dd10<- data.table(DD1040_01 =ifelse(Tave3[[ 1 ]] <= 18 , 302.126 /(1 + exp(-(Tave3[[ 1 ]]- 14.71 )/ 2.3879 )), 32.1984 *Tave3[[ 1 ]]+ -340 ));head(dd10)
      dd10$DD1040_02 =ifelse(Tave3[[ 2 ]] <= 18 , 270.2587 /(1 + exp(-(Tave3[[ 2 ]]- 14.78 )/ 2.4094 )), 26.4552 *Tave3[[ 2 ]]+ -250 )
      dd10$DD1040_03 =ifelse(Tave3[[ 3 ]] <= 16 , 229.4616 /(1 + exp(-(Tave3[[ 3 ]]- 13.46 )/ 2.0635 )), 30.2728 *Tave3[[ 3 ]]+ -300 )
      dd10$DD1040_04 =ifelse(Tave3[[ 4 ]] <= 16 , 387.3743 /(1 + exp(-(Tave3[[ 4 ]]- 16.39 )/ 2.9742 )), 29.7395 *Tave3[[ 4 ]]+ -300 )
      dd10$DD1040_05 =ifelse(Tave3[[ 5 ]] <= 16 , 251.5565 /(1 + exp(-(Tave3[[ 5 ]]- 13.79 )/ 2.6216 )), 31.231 *Tave3[[ 5 ]]+ -320 )
      dd10$DD1040_06 =ifelse(Tave3[[ 6 ]] <= 13.7 , 183.1119 /(1 + exp(-(Tave3[[ 6 ]]- 12.74 )/ 1.7739 )), 28.8581 *Tave3[[ 6 ]]+ -280 )
      dd10$DD1040_07 =ifelse(Tave3[[ 7 ]] <= 15 , 202.1216 /(1 + exp(-(Tave3[[ 7 ]]- 13.13 )/ 1.72 )), 30.7995 *Tave3[[ 7 ]]+ -310 )
      dd10$DD1040_08 =ifelse(Tave3[[ 8 ]] <= 14 , 178.9318 /(1 + exp(-(Tave3[[ 8 ]]- 12.65 )/ 1.6275 )), 30.3287 *Tave3[[ 8 ]]+ -300 )
      dd10$DD1040_09 =ifelse(Tave3[[ 9 ]] <= 14 , 183.2781 /(1 + exp(-(Tave3[[ 9 ]]- 12.54 )/ 2.2868 )), 29.4147 *Tave3[[ 9 ]]+ -290 )
      dd10$DD1040_10 =ifelse(Tave3[[ 10 ]] <= 14 , 487.0833 /(1 + exp(-(Tave3[[ 10 ]]- 17.43 )/ 3.2032 )), 30.2894 *Tave3[[ 10 ]]+ -300 )
      dd10$DD1040_11 =ifelse(Tave3[[ 11 ]] <= 15 , 229.3626 /(1 + exp(-(Tave3[[ 11 ]]- 13.69 )/ 2.244 )), 29.6914 *Tave3[[ 11 ]]+ -300 )
      dd10$DD1040_12 =ifelse(Tave3[[ 12 ]] <= 14 , 159.4986 /(1 + exp(-(Tave3[[ 12 ]]- 12.21 )/ 1.7726 )), 29.8672 *Tave3[[ 12 ]]+ -300 )
      dd10$DD1040_12py =ifelse(Tave12 <= 14 , 159.4986 /(1 + exp(-(Tave12[[1]]- 12.21 )/ 1.7726 )), 29.8672 *Tave12[[1]]+ -300 );head(dd10)
      
      dd10[dd10<0] <- 0; head(dd10)
      dd10$DD1040_wt <- rowSums(dd10[,c(1,2,13)])
      dd10$DD1040_sp <- rowSums(dd10[,3:5])
      dd10$DD1040_sm <- rowSums(dd10[,6:8])
      dd10$DD1040_at <- rowSums(dd10[,9:11])
      dd10$DD1040 <- rowSums(dd10[,1:12])
      dd10 <- round(dd10,0);head(dd10)
      dd10a <- data.table(dd10[,c(1:12,14:18)]);dd10a
      DT0_ecr[row0:row2,'DD1040' := dd10a$DD1040]; DT0_ecr
      rm(dd10,dd10a);gc()
      
      ##CMI ---------------------------------------------------
      El2 = xyl_set$elev0; El2 =ifelse(El2== -9999, xyl_set$elev, El2)
      i=1
      for(i in 1: 12){
        e_tmax_kpa = 0.61078 * exp(17.269 * Tmax3[[i]] / (237.3 + Tmax3[[i]]))
        e_tmin_kpa = 0.61078 * exp(17.269 * Tmin3[[i]] / (237.3 + Tmin3[[i]]))
        e_tdew_kpa = 0.61078 * exp(17.269 * (Tmin3[[i]] - 2.5) / (237.3 + Tmin3[[i]] - 2.5))
        vpd_kpa = 0.5 * (e_tmax_kpa + e_tmin_kpa) - e_tdew_kpa
        tmean_c = 0.5 * (Tmax3[[i]] + Tmin3[[i]])
        k_trf = (tmean_c + 5) / 15
        k_trf[k_trf > 1] <- 1
        k_trf[k_trf < 0] <- 0
        pet_mm = 93 * vpd_kpa * k_trf * exp(El2 / 9300)
        cmi0 = data.table((ppt3[[i]] - pet_mm) / 10)
        if(i==1){CMI=cmi0}else{CMI=data.table(CMI,cmi0)}
        #if(i==1){pet=pet_mm}else{pet=data.table(pet,pet_mm)}
      }
      
      #prev 12
      e_tmax_kpa = 0.61078 * exp(17.269 * Tmax12 / (237.3 + Tmax12))
      e_tmin_kpa = 0.61078 * exp(17.269 * Tmin12 / (237.3 + Tmin12))
      e_tdew_kpa = 0.61078 * exp(17.269 * (Tmin12 - 2.5) / (237.3 + Tmin12 - 2.5))
      vpd_kpa = 0.5 * (e_tmax_kpa + e_tmin_kpa) - e_tdew_kpa
      tmean_c = 0.5 * (Tmax12 + Tmin12)
      k_trf = (tmean_c + 5) / 15
      k_trf[k_trf > 1] <- 1
      k_trf[k_trf < 0] <- 0
      pet_mm = 93 * vpd_kpa * k_trf * exp(El2 / 9300)
      cmi12py = (ppt3[[i]] - pet_mm) / 10
      CMI$CMI_12py <- cmi12py
      
      names(CMI)[1:12] <- paste0('CMI',substr(101:112,2,3));head(CMI)
      #names(pet)[1:12] <- paste0('pet',substr(101:112,2,3));head(pet)
      CMI$CMI_wt <- rowSums(CMI[,c(1,2,13)]);
      CMI$CMI_sp <- rowSums(CMI[,3:5]);
      CMI$CMI_sm <- rowSums(CMI[,6:8]);
      CMI$CMI_at <- rowSums(CMI[,9:11]);
      CMI$CMI <- rowSums(CMI[,1:12]);
      CMI <- round(CMI*100,0);head(CMI)
      CMI2 <- data.table(CMI[,c(1:12,14:18)]);CMI2
      DT2_cmi[row0:row2,names(CMI2) := CMI2];DT2_cmi
      #DT2_cmi[row0:row2,names(pet) := round(pet*10,2)];DT2_cmi
      rm(CMI,CMI2)
      rm(yVar,nffd);gc()
      
      # rsds --------------------------------
      # do nothing for Ref or GCM --
      if(grepl('His',prd_class)){
        # message('gereate rsds')
        td03 = (DT1_ttr[row0:row2]$Tmax03 - DT1_ttr[row0:row2]$Tmin03)/10;head(td03)
        rsds <- data.table(rsds01 = 29.99 + (-0.716 * xyl_set$y) + (-0.185 * DT0_ecr[row0:row2]$ErefR01) + (0.00326 * xyl_set$elev) + (0.00406 * xyl_set$y^2) + (0.000787 * DT0_ecr[row0:row2]$ErefR01^2) + (-0.000000348 * xyl_set$elev^2) + (0.00357 * xyl_set$y * DT0_ecr[row0:row2]$ErefR01) + (-0.0000480 * xyl_set$y * xyl_set$elev) + (-0.0000278 * DT0_ecr[row0:row2]$ErefR01 * xyl_set$elev) + (0.00000150 * xyl_set$y * DT0_ecr[row0:row2]$ErefR01 * xyl_set$elev))
        rsds$rsds02 = -8.66 + (0.206 * xyl_set$y) + (-0.491 * xyl_set$x) + (0.473 * DT0_ecr[row0:row2]$RH02) + (0.000856 * xyl_set$y^2) + (0.000169 * xyl_set$x^2) + (-0.000741 * DT0_ecr[row0:row2]$RH02^2) + (0.00893 * xyl_set$y * xyl_set$x) + (-0.00718 * xyl_set$y * DT0_ecr[row0:row2]$RH02) + (0.00703 * xyl_set$x * DT0_ecr[row0:row2]$RH02) + (-0.000119 * xyl_set$y * xyl_set$x * DT0_ecr[row0:row2]$RH02)
        rsds$rsds03 = 53.4 + (-0.613 * xyl_set$y) + (-3.72 * td03) + (0.432 * xyl_set$x) + (-0.00168 * xyl_set$y^2) + (-0.0112 * td03^2) + (0.000203 * xyl_set$x^2) + (0.0680 * xyl_set$y * td03) + (-0.00620 * xyl_set$y * xyl_set$x) + (-0.0472 * td03 * xyl_set$x) + (0.000769 * xyl_set$y * td03 * xyl_set$x)
        rsds$rsds04 = -57.2 + (1.18 * xyl_set$y) + (1.18 * DT0_ecr[row0:row2]$RH04) + (-1.00 * xyl_set$x) + (-0.000156 * xyl_set$y^2) + (-0.000715 * DT0_ecr[row0:row2]$RH04^2) + (0.000147 * xyl_set$x^2) + (-0.0179 * xyl_set$y * DT0_ecr[row0:row2]$RH04) + (0.0160 * xyl_set$y * xyl_set$x) + (0.0142 * DT0_ecr[row0:row2]$RH04 * xyl_set$x) + (-0.000220 * xyl_set$y * DT0_ecr[row0:row2]$RH04 * xyl_set$x)
        rsds$rsds05 = 21.9 + (-0.0742 * DT0_ecr[row0:row2]$ErefR05) + (0.00554 * xyl_set$elev) + (-0.00834 * xyl_set$x) + (0.000243 * DT0_ecr[row0:row2]$ErefR05^2) + (0.000000241 * xyl_set$elev^2) + (-0.000211 * xyl_set$x^2) + (-0.000138 * DT0_ecr[row0:row2]$ErefR05 * xyl_set$elev) + (-0.000472 * DT0_ecr[row0:row2]$ErefR05 * xyl_set$x) + (0.0000479 * xyl_set$elev * xyl_set$x) + (-0.00000135 * DT0_ecr[row0:row2]$ErefR05 * xyl_set$elev * xyl_set$x)
        rsds$rsds06 = 83.5 + (-1.32 * DT0_ecr[row0:row2]$RH06) + (-0.0587 * xyl_set$x) + (-1.48 * xyl_set$y) + (0.00314 * DT0_ecr[row0:row2]$RH06^2) + (-0.000728 * xyl_set$x^2) + (-0.00205 * xyl_set$y^2) + (-0.00183 * DT0_ecr[row0:row2]$RH06 * xyl_set$x) + (0.0240 * DT0_ecr[row0:row2]$RH06 * xyl_set$y) + (-0.00745 * xyl_set$x * xyl_set$y) + (0.000110 * DT0_ecr[row0:row2]$RH06 * xyl_set$x * xyl_set$y)
        rsds$rsds07 = 7.70 + (-0.00986 * DT0_ecr[row0:row2]$ErefR07) + (-0.198 * xyl_set$x) + (0.0158 * xyl_set$elev) + (0.000164 * DT0_ecr[row0:row2]$ErefR07^2) + (-0.000995 * xyl_set$x^2) + (0.000000309 * xyl_set$elev^2) + (-0.000159 * DT0_ecr[row0:row2]$ErefR07 * xyl_set$x) + (-0.000190 * DT0_ecr[row0:row2]$ErefR07 * xyl_set$elev) + (0.000143 * xyl_set$x * xyl_set$elev) + (-0.00000177 * DT0_ecr[row0:row2]$ErefR07 * xyl_set$x * xyl_set$elev)
        rsds$rsds08 = 5.79 + (0.0532 * DT0_ecr[row0:row2]$ErefR08) + (-0.0923 * xyl_set$x) + (0.00843 * xyl_set$elev) + (-0.000157 * DT0_ecr[row0:row2]$ErefR08^2) + (-0.000654 * xyl_set$x^2) + (0.000000281 * xyl_set$elev^2) + (-0.000533 * DT0_ecr[row0:row2]$ErefR08 * xyl_set$x) + (-0.000138 * DT0_ecr[row0:row2]$ErefR08 * xyl_set$elev) + (0.0000691 * xyl_set$x * xyl_set$elev) + (-0.00000126 * DT0_ecr[row0:row2]$ErefR08 * xyl_set$x * xyl_set$elev)
        rsds$rsds09 = 4.90 + (0.102 * DT0_ecr[row0:row2]$ErefR09) + (-0.000119 * xyl_set$elev) + (-0.00471 * xyl_set$x) + (-0.000457 * DT0_ecr[row0:row2]$ErefR09^2) + (0.000000371 * xyl_set$elev^2) + (-0.000133 * xyl_set$x^2) + (-0.0000745 * DT0_ecr[row0:row2]$ErefR09 * xyl_set$elev) + (-0.000686 * DT0_ecr[row0:row2]$ErefR09 * xyl_set$x) + (-0.00000248 * xyl_set$elev * xyl_set$x) + (-0.000000721 * DT0_ecr[row0:row2]$ErefR09 * xyl_set$elev * xyl_set$x)
        rsds$rsds10 = -0.191 + (0.0901 * DT0_ecr[row0:row2]$ErefR10) + (0.216 * xyl_set$y) + (0.000974 * xyl_set$elev) + (-0.00262 * xyl_set$y^2) + (0.000000234 * xyl_set$elev^2) + (0.000349 * DT0_ecr[row0:row2]$ErefR10 * xyl_set$y) + (-0.0000175 * DT0_ecr[row0:row2]$ErefR10 * xyl_set$elev) + (-0.0000160 * xyl_set$y * xyl_set$elev) + (0.000000964 * DT0_ecr[row0:row2]$ErefR10 * xyl_set$y * xyl_set$elev)
        rsds$rsds11= 11.6 + (-0.184 * xyl_set$y) + (-0.0492 * DT0_ecr[row0:row2]$ErefR11) + (0.00265 * xyl_set$elev) + (0.000313 * xyl_set$y^2) + (0.000468 * DT0_ecr[row0:row2]$ErefR11^2) + (-0.000000115 * xyl_set$elev^2) + (0.00272 * xyl_set$y * DT0_ecr[row0:row2]$ErefR11) + (-0.0000391 * xyl_set$y * xyl_set$elev) + (-0.0000282 * DT0_ecr[row0:row2]$ErefR11 * xyl_set$elev) + (0.00000144 * xyl_set$y * DT0_ecr[row0:row2]$ErefR11 * xyl_set$elev)
        rsds$rsds12 = 20.5 + (-0.482 * xyl_set$y) + (-0.110 * DT0_ecr[row0:row2]$ErefR12) + (0.00364 * xyl_set$elev) + (0.00259 * xyl_set$y^2) + (0.000584 * DT0_ecr[row0:row2]$ErefR12^2) + (-0.000000274 * xyl_set$elev^2) + (0.00316 * xyl_set$y * DT0_ecr[row0:row2]$ErefR12) + (-0.0000542 * xyl_set$y * xyl_set$elev) + (-0.0000293 * DT0_ecr[row0:row2]$ErefR12 * xyl_set$elev) + (0.00000130 * xyl_set$y * DT0_ecr[row0:row2]$ErefR12 * xyl_set$elev)
        rsds$rsds12py = rsds$rsds12
        #previous year ---
        rsds$rsds12py = 20.5 + (-0.482 * xyl_set$y) + (-0.110 * ErefR12py) + (0.00364 * xyl_set$elev) + (0.00259 * xyl_set$y^2) + (0.000584 * ErefR12py^2) + (-0.000000274 * xyl_set$elev^2) + (0.00316 * xyl_set$y * ErefR12py) + (-0.0000542 * xyl_set$y * xyl_set$elev) + (-0.0000293 * ErefR12py * xyl_set$elev) + (0.00000130 * xyl_set$y * ErefR12py * xyl_set$elev)
        rsds[rsds<0] <- 0; head(rsds)
        rsds$rsds_wt = rowMeans(rsds[,c(13,1,2)])
        rsds$rsds_sp = rowMeans(rsds[,c(3:5)])
        rsds$rsds_sm = rowMeans(rsds[,c(6:8)])
        rsds$rsds_at = rowMeans(rsds[,c(9:11)])
        rsds$rsds = rowMeans(rsds[,c(1:12)])
        rsds <- round(rsds,1);rsds
        rsds2 <- data.table(rsds)
        DT1_ttr[row0:row2,names(rsds2) := rsds2*10];DT1_ttr
        rm(td03,rsds,rsds2);
        rm(e_tdew_kpa,e_tmax_kpa,e_tmin_kpa,Tave12,Tave3,Tmax_h,Tmax3,Tmax12);gc()
        rm(Tmin12,Tmin3,Tminx,PPT12,ppt_h,ppt3);gc()
      }
      message(paste0(floor(Set/nSet*100), '% completed for ',prd))
    } #end of subset looping
    ##combining=========
    DT <- data.table(DT0_pdd,DT0_ddnp,DT0_ecr,DT1_ttr/10,DT2_cmi/100);head(DT)
    
    ## outPut ============================================================================>
    message('Saving files ...')
    VL=1
    if(length(varList)==1){VL=ifelse(varList=='Y'|varList=='M'|varList=='S'|varList=='PM'|varList=='DM'|varList=='YM'|varList=='MY'|varList=='YS'|varList=='SY'|varList=='YSM'|varList=='MSY'|varList=='YMS',0,1);VL}
    if(VL==0){if(varList=='Y'){varList=varList_Y;VL=1}}
    if(VL==0){if(varList=='S'){varList=varList_S;VL=1}}
    if(VL==0){if(varList=='PM'){varList=varList_PM;VL=1}}
    if(VL==0){if(varList=='DM'){varList=varList_DM;VL=1}}
    if(VL==0){if(varList=='M'){varList=c(varList_PM,varList_DM);VL=1}}
    if(VL==0){if(varList=='YS'|varList=='SY'){varList=c(varList_Y,varList_S);VL=1}}
    if(VL==0){if(varList=='YM'|varList=='MY'){varList=c(varList_Y,varList_PM,varList_DM);VL=1}}
    if(VL==0){if(varList=='YSM'|varList=='MSY'|varList=='YMS'){varList=c(varList_Y,varList_S,varList_PM,varList_DM)}}
    
    #output file path/names -------------
    if(is.character(inputFile)){
      input_splt <- strsplit(inputFile,'/')[[1]];input_splt
      inF_name = input_splt[length(input_splt)];inF_name
      inF_name2 = substr(inF_name,1,(nchar(inF_name)-4));inF_name2
      prd_name = gsub('/','_',substr(prd,1,(nchar(prd)-4)));prd_name
    }
    if(inputFile_class=="CSV"){
      out_csv <- data.table(input,DT[,..varList]);head(out_csv)
      outF <- paste0(outDir,inF_name2,'_',gsub('/','_',substr(prd,1,(nchar(prd)-4))),'.csv');outF
      data.table::fwrite(out_csv,outF)
    }# end of CSV output
    
    if(inputFile_class=="dataFrame"){
      out_csv <- data.table(input,DT[,..varList]);head(out_csv)
      inF_name2 <- prefix
      outF <- paste0(outDir,inF_name2,'_',gsub('/','_',substr(prd,1,(nchar(prd)-4))),'.csv');outF
      print(outF)
      data.table::fwrite(out_csv,outF)
    }# end of CSV output
    
    if(inputFile_class=="Raster"){
      # # Create output directory structure
      out_subDir <- file.path(outDir, inF_name2)
      out_subDir2 <- file.path(out_subDir, prd_name)
      dir.create(out_subDir, showWarnings = FALSE)
      dir.create(out_subDir2, showWarnings = FALSE)
      
      # Process one variable at a time
      # i=1
      for(i in 1:length(varList)) {
        # Create raster using coordinates and values
        temp_df <- data.frame(x=xyl_input[,1], y=xyl_input[,2], value=DT[[varList[i]]])
        temp_rast <- terra::rast(temp_df)
        names(temp_rast) <- varList[i]
        # Set extent and origin
        terra::ext(temp_rast) <- ext_input
        terra::origin(temp_rast) <- origin_input
        # Write to file
        outF <- file.path(out_subDir2, paste0(varList[i], '.tif'))
        writeRaster(temp_rast, outF, overwrite=TRUE)
        # Optional: check the output
        # plot(temp_rast, main=varList[i])
        # Clean up
        rm(temp_rast, temp_df)
        gc()
      }
      rm(DT);    gc()
    }#end of tif output
    #rm(stk_out);gc()
    message(prd, ' --- completed')
  } #end of periodList looping
  
  ## Return--------
  #if(inputFile_class=="CSV"|inputFile_class=="dataFrame"){return(out_csv)}
  #if(inputFile_class=="Raster"){return(stk_out)}
  print(paste0('Completed for ', length(periodList),' periods'))
  rm(list=ls());gc()
} # end of ClimateNAr <===========================================================================