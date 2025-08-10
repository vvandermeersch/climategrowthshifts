This README file was generated on 2025-07-18 by Pieter Zuidema.
GENERAL INFORMATION

1. Title of Dataset: Pantropical tree rings show small effects of drought on stem growth
2. Author Information
   A. Principal Investigator Contact Information
   Name: Pieter Zuidema
   Institution: Wageningen University
   Address: Wageningen, the Netherlands
   Email: [pieter.zuidema@wur.nl](mailto:pieter.zuidema@wur.nl)
   B. Associate or Co-investigator Contact Information
   Name: Flurin Babst
   Institution: University of Arizona
   Address: Tucson, AZ, USA
   Email: [babst@arizona.edu](mailto:babst@arizona.edu)
3. Date of data collection: compiled tree-ring data from different sources. Tree-ring data
4. Geographic location of data collection: pantropical (30 degree N -30 degree S)
   SHARING/ACCESS INFORMATION
5. Licenses/restrictions placed on the data: CC0 1.0 Universal (CC0 1.0) Public Domain
6. Links to publication that cite or use the data:
   Zuidema, P.A. et al. (2025). Pantropical tree rings show small effects of drought on stem growth. Science, DOI: 10.1126/science.adq6607
7. Links to other publicly accessible locations of the data: Raw tree-ring width measurements for all 483 chronologies are publicly available free of charge and for all research purposes here:
   National Centers for Environmental Information, National Oceanic and Atmospheric Administration, International tree- ring data bank (ITRDB) (NOAA, 2025); [https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring](https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring)
8. Links/relationships to ancillary data sets: None
9. Was data derived from another source? Yes.
   A. If yes, list source(s): Raw tree-ring width data for 339 chronologies were obtained from the ITRDB (International Tree-Ring Databank, [https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring](https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring)); those for 144 chronologies from co-authors who contributed to the Tropical Tree-ring Network network (see [www.tropicaltreeringnetwork.org](http://www.tropicaltreeringnetwork.org)).
10. Recommended citation for this dataset:
    Zuidema, P.A. et al. (2025). Data from: Pantropical tree rings show small effects of drought on stem growth. Dryad. [https://doi.org/10.5061/dryad.hx3ffbgq4](https://doi.org/10.5061/dryad.hx3ffbgq4)
    DATA & FILE OVERVIEW
11. File List included in Data.zip:
    In /Data folder:
    A) metadata_FINAL.csv
    B) Relative_climate_densities.csv
    C) plots_meta.csv
    D) plots_additional_mortality.csv
    E) sea_1_standard.Rdata
    F) sea_2_wet.csv
    G) sea_3_series.csv
    H) sea_4_splines.csv
    I) sea_5_CRU_ERA.Rdata,
    J) sea_6_CRU_TC.Rdata
    K) sea_7_fivepercent.csv
    L) sea_8_repeated.csv
    M) MODIS_tree_cover.tif
    N) S_mean_raster.tif
    In /Data/1_Climate folder:
    O) 6 txt files for 3 climate variables (P, VPD, CWD) and two seasons (dry, wet)
    P) 6 txt files for 3 climate variables (P, VPD, CWD) and two seasons (dry, wet)
    Q) 6 txt files for 3 climate variables (P, VPD, CWD) and two seasons (dry, wet)
    R) 5 txt files with chronologies for 5 detrending methods.
    In /Data/Interpolation folder:
    S) 24 GeoTiff files with climatic interpolation output
    T) 4  GeoTiff files with maps based on interpolation output
    In /Data/Series_Level folder
    U) 483 netcdf files with names corresponding to SiteID in the metadata_FINAL.csv file
    In /Data/Artifical_Growth_Reduction folder
    V) Randomyear.txt
    W) 3 txt files with articially perturbed chronologies
12. Relationship between files, if important: None
13. Additional related data collected that was not included in the current data package: Yes, (A) raw tree-ring width data for all 483 chronologies is not included in the datapackage but can be retrieved from the ITRDB (International Tree-Ring Databank, [https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring](https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring)); (B) raw climate data is not included in this data package; this was collected and is available from four different products (WorldClim, CRU, ERA5, TerraClimate).
14. Are there multiple versions of the dataset? No
    A. If yes, name of file(s) that was updated: NA
    i. Why was the file updated? NA
    ii. When was the file updated? NA

    **DATA-SPECIFIC INFORMATION FOR: metadata_FINAL.csv**
15. Description: Metadata on chronologies: location, species, length, climate.
16. Number of variables: 34
17. Number of cases/rows: 483
18. Variable List:

* SiteID: Unique identifier of the chronology
* SiteName: Name of the site where trees for chronology were sampled
* Source: Whether data was obtained from co-authors for this study and now added to the ITRDB (ITRDB_added_in_2025), had been previously obtained from contributors to the network (ITRDB_added_in_2022), or had already been in the ITRDB and was obtained from there (ITRDB_before_2022)
* Country: Country where chronology trees for the chronology were sampled
* Region: One of three tropical regions
* Sp_code: Abbreviation of the sampled species name consisting of first 2 letters of genus and 3 first letters of the species epithet
* Genus: genus name
* Sp: species epithet name
* lat: Latitude of sampling site (degrees)
* lon: Longitude of sampling site (degrees)
* Elevation: Elevation of sampling site (m above sea level)
* latCRU: Latitude of site used to assign site to a CRU climate pixel (degrees)
* lonCRU: Longitude of site used to assign site to a CRU climate pixel (degrees)
* Shulman: Whether ring dating for Southern Hemisphere follows the Shulman convention, which stipulates that the ring is dated as the year when ring formulation started. Yes: follows convention, No: does not follow convetion, NA: not applicable (for Northern Hemisphere sites)
* ntrees: number of trees included in the chronology
* nseries: number of series included in the chronology
* Start: start year of the chronology
* End: end year of the chronology
* Rbar:
* MATwclim: mean annual temperature for the site from WorldClim (degree C)
* MAPwclim: mean annual precipitation for the site from WorldClim (mm/y)
* Wet_start: month number (Jan=1; Dec=12) when the wet season starts at the site
* Wet_end: month number (Jan=1; Dec=12) when the wet season ends at the site
* Dry_start: month number (Jan=1; Dec=12) when the dry season starts at the site
* Dry_end: month number (Jan=1; Dec=12) when the dry season ends at the site
* Nseasons: number of seasons at the site
* Wetmonths: number of wet months at the site
* GymAngio: whether study species is a Gymnosperm or Angiosperm
* Family: taxonomic family of the sampled species
* WD: wood density of the sampled species (kg/m3)
* studyCode: unique identifier of the tree-ring dataset in the ITRDB
* doi: unique link to the tree-ring dataset in the ITRDB
* NOAAStudyId: unique identifier of the tree-ring dataset in the ITRDB
* NOAASiteId: unique identifier of the study site in the ITRDB

1. Missing data codes: NA
2. Specialized formats or other abbreviations used: None

   **DATA-SPECIFIC INFORMATION FOR: Relative_climate_densities.csv**
3. Description: Relative climate densities of chronologies for analyses weighted for climate representativeness. Used to produce Figure S3
4. Number of variables: 4
5. Number of cases/rows: 483
6. Variable List:

* SiteID: Unique identifier of the chronology
* Rel_density_MAP: The relative density of tropical land area with woody vegetation (tree cover >10%) for the mean annual precipitation (MAP) at the chronology site. Ranges from 0 (no tropical land area with MAP equal to site) to 1 (MAP at site is within most abundant MAP category for tropical land area).
* Rel_density_CWD: The relative density of tropical land area with woody vegetation (tree cover >10%) for the climatic water deficit (CWD) at the chronology site. Ranges from 0 (no tropical land area with CWD equal to site) to 1 (CWD at site is within most abundant CWD category for tropical land area).
* Rel_density_VPD: The relative density of tropical land area with woody vegetation (tree cover >10%) for the vapor pressure deficit (VPD) at the chronology site. Ranges from 0 (no tropical land area with VPD equal to site) to 1 (VPD at site is within most abundant VPD category for tropical land area).

1. Missing data codes: None
2. Specialized formats or other abbreviations used: None

   **DATA-SPECIFIC INFORMATION FOR: plots_meta.csv**
3. Description: Metadata and coefficients of growth-mortality associations for permanent sample plots and coefficients used to estimate drought-induced mortality. Data obtained from Russo et al 2021, ref 49 in paper. Used to produce Figure S18.
4. Number of variables: 4
5. Number of cases/rows: 10
6. Variable List:

* SiteID: Unique identifier of the plot
* MAPwclim: mean annual precipitation for the plot from WorldClim (mm/y)
* MATwclim: mean annual temperature for the plot from WorldClim (degree C)
* Beta_0: the intercept of the growth-mortality relationship in Russo et al 2021 for the plot
* Beta_2: the coefficient defining the growth-mortality relationship in Russo et al. 2021 for the plot
* R2: R2 value of the growth-mortality relationship in Russo et al. 2021 for the plot
* Nspp: number of species on which the average growth-mortality relationship is based
* SiteID_R2_nspp: a text compilation of the SiteID, R2 and nspp to be shown in Figure S18.

1. Missing data codes: None
2. Specialized formats or other abbreviations used: None

   **DATA-SPECIFIC INFORMATION FOR: plots_additional_mortality.csv**
3. Description: Metadata and coefficients of growth-mortality associations for permanent sample plots and coefficients used to estimate drought-induced mortality. Data obtained from Russo et al 2021, ref 49 in paper. Used to produce Figure S18.
4. Number of variables: 4
5. Number of cases/rows: 10
6. Variable List:

* SiteID: Unique identifier of the plot. ALL for all plots together.
* Beta0_of_site: the intercept of the growth-mortality relationship in Russo et al 2021 for the plot
* Beta2_of_site: the coefficient defining the growth-mortality relationship in Russo et al. 2021 for the plot
* Beta1_adjusted_per_site: coefficient of the size-mortality relationship, which was not provided in Russo et al 2021 and therefore estimated to obtain a 1-% annual tree mortality across all sites.
* Dij: diameter of the tree, set to 200 mm
* DBH_growth: normal stem growth in diameter at breast height, set to 2 mm/y
* Tauij: power transformation of growth rate, calculated as DBH_growth^(0.45)
* Logit_mort_orig: logit of 5-year mortality rate for normal growth
* Mort_py_orig: annual mortality rate for normal growth
* Grwth_reduct: fraction of growth reduction applied to calculate new mortality rate
* Logit_mort_NEW: logit of 5-year mortality rate when growth is reduced by factor Grwth_reduct
* Mort_py_NEW: annual mortality rate when growth is reduced by factor Grwth_reduct
* Mort_increase: increase in mortality rate due to growth reduction by Grwth_reduct
* R2_site: R2 value of the growth-mortality relationship in Russo et al. 2021 for the plot
* n_spp: number of species on which the average growth-mortality relationship is based
* SiteID_R2_nspp: a text compilation of the SiteID, R2 and nspp to be shown in Figure S18.
* Mort_increase_minSD: increase in mortality rate due to growth reduction by Grwth_reduct when Beta2 was reduced by one standard deviation (SD=0.802) to account for across-site variation.
* Mort_increase_plusSD: increase in mortality rate due to growth reduction by Grwth_reduct when Beta2 was increased by one standard deviation (SD=0.802) to account for across-site variation.

1. Missing data codes: None
2. Specialized formats or other abbreviations used: None

   **DATA-SPECIFIC INFORMATION FOR: sea_1_standard.Rdata**
3. Description: Results of Superposed Epoch Analysis for 10% driest years. Used to produce several main and SI Figures.
4. Number of variables: 28
5. Number of cases/rows: 8166
6. Variable List:

* lag: lag of drought year (0=drought year)
* se: the scaled mean RWI of the drought years
* se.unscaled: the unscaled mean RWI of the drought years
* p: p-value of SEA analysis
* ci.95.lower: lower limit of the 95% confidence interval for the scaled RWI of drought years
* ci.95.upper: upper limit of the 95% confidence interval for the scaled RWI of drought years
* SiteID: Unique identifier of the chronology
* clim_var: drought type (P=preciptation, VPD=vapor pressure deficit, CWD=climatic water deficit)
* season: droughts during the dry or wet season
* spline: detrending flexibility, the period (y) at which 50% frequency cut-off is applied
* nr_evnts: number of drought events in the chronology
* evnts: list of calendar years of drought events
* evnts_clim: climate (P, VPD or CWD) during specified season (wet, dry) of dry years
* ar1: first-order autocorrelation in chronology
* anom_chr: list of RWI anomalies during drought years
* mn_rwi: mean ring width index of normal years
* mn_clim: mean climate (P, VPD or CWD) during specified season (wet, dry) of normal years
* sd_clim: standard deviation of climate (P, VPD or CWD) during specified season (wet, dry) of normal years
* iav_clim: interannual variability of climate (P, VPD or CWD) during specified season (wet, dry) of normal years
* mn_evnts_clim: mean climate (P, VPD or CWD) during specified season (wet, dry) of dry years
* mn_lag1_clim: mean climate (P, VPD or CWD) during specified season (wet, dry) of years with lag=1
* mn_lag2_clim: mean climate (P, VPD or CWD) during specified season (wet, dry) of years with lag=2
* sign: whether SEA is significant (1: p<0.05, 0: p>0.05)
* RWI_anomaly: mean RWI anomaly during drought years, calculated as se.unscaled - mn_rwi
* times_sd_clim: climate anomaly during drought years as the number of standard deviations (SD) away from the mean climate for the specified drought type (P, VPD, CWD) during the specified season (wet, dry)
* times_sd_clim_lag1: climate anomaly during lag-1 years as the number of standard deviations (SD) away from the mean climate for the specified drought type (P, VPD, CWD) during the specified season (wet, dry)
* times_sd_clim_lag2: climate anomaly during lag-2 years as the number of standard deviations (SD) away from the mean climate for the specified drought type (P, VPD, CWD) during the specified season (wet, dry)
* clim_var_long: drought type, full names

  **DATA-SPECIFIC INFORMATION FOR: sea_2_wet.csv**

1. Description: Results of Superposed Epoch Analysis for 10% wettest years. Used to produce main Figure 3 and SI Figure S7B & S13.
2. Number of variables: 16
3. Number of cases/rows: 8253
4. Variable List:

* subset of variables sea_1_standard.Rdata, but then for wet extreme years instead of drought years

  **DATA-SPECIFIC INFORMATION FOR: sea_3_series.csv**

1. Description: Results of Superposed Epoch Analysis for ring-width series. Used to produce SI Figure S8.
2. Number of variables: 13
3. Number of cases/rows: 7758
4. Variable List:

* propsign: proportion of series for which SEA was significant
* propsign_neg: proportion of series for which SEA was negative and significant
* propsign_pos: proportion of series for which SEA was positive and significant
* SiteID: Unique identifier of the chronology
* clim_var: drought type (P=preciptation, VPD=vapor pressure deficit, CWD=climatic water deficit)
* season: droughts during the dry or wet season
* spline: detrending flexibility, the period (y) at which 50% frequency cut-off is applied
* nr_evnts: number of drought events in the chronology
* N_series: number of series in the chronology
* N_series_analysed: number of series in the chronology for which SEA analysis could be done
* clim_var_long: drought type, full names

  **DATA-SPECIFIC INFORMATION FOR: sea_4_splines.csv**

1. Description: Results of Superposed Epoch Analysis with 5 detrending procedures. Used to produce SI Figure S9
2. Number of variables: 16
3. Number of cases/rows: 33012
4. Variable List:

* subset of variables from sea_1_standard.Rdata

  **DATA-SPECIFIC INFORMATION FOR: sea_5_CRU_ERA.Rdata**

1. Description: Results of Superposed Epoch Analysis using an alternative climate product (ERA5) in addition to CRU. Used to produce SI Figures S10 and S11.
2. Number of variables: 16
3. Number of cases/rows: 16506
4. Variable List:

* subset of variables from sea_1_standard.Rdata
* clim_product: which climate product was used for the SEA analysis: CRU or ERA5

  **DATA-SPECIFIC INFORMATION FOR: sea_6_CRU_TC.Rdata**

1. Description: Results of Superposed Epoch Analysis using an alternative climate product (TerraClimate, TC) in addition to CRU. Used to produce SI Figures S10 and S11.
2. Number of variables: 16
3. Number of cases/rows: 15744
4. Variable List:

* subset of variables from sea_1_standard.Rdata
* clim_product: which climate product was used for the SEA analysis: CRU or TerraClimate (TC)

  **DATA-SPECIFIC INFORMATION FOR: sea_7_fivepercent.csv**

1. Description: Results of Superposed Epoch Analysis for the 5% driest years. Used to produce SI Figure S12A.
2. Number of variables: 16
3. Number of cases/rows: 7362
4. Variable List:

* subset of variables from sea_1_standard.Rdata
  \
  **DATA-SPECIFIC INFORMATION FOR: sea_8_repeated.csv**

1. Description: Results of Superposed Epoch Analysis for repeated (2-year) droughts. Used to produce SI Figure S12C.
2. Number of variables: 16
3. Number of cases/rows: 8253
4. Variable List:

* subset of variables from sea_1_standard.Rdata

  DATA-SPECIFIC INFORMATION FOR: MODIS_tree_cover.tif

1. Description: Global treecover SpatRaster map. Used to produce background tree cover in all maps.
2. Number of variables: 1
3. Number of rows, columns: 900, 2160
4. Variable List:

* Tree cover proportion, only when >0.1

  **DATA-SPECIFIC INFORMATION FOR: S_mean_raster.tif**

1. Description: Global SpatRaster map of tree species richness from Liang et al. (2022). Used to produce Fig S1E. C.
2. Number of variables: 1
3. Number of rows, columns: 7200, 14400
4. Variable List:

* Estimated number of tree species (/ha)

  **DATA-SPECIFIC INFORMATION FOR: 6 files with names CRU__**_100mm_TROTS4.txt**

1. Description: CRU climate data for all 483 sites. Used as input for all SEA analyses. File names indicate climate variable (pre=P, vpd=VPD, cwd=CWD).
2. Number of variables: 483
3. Number of rows: 120
4. Variable List:

* Each column contains the seasonal climate data for one chronology (site). Column name corresponds to the SiteID from metadata_FINAL.csv
* Each row corresponds to a year, with the row name representing the calendar year

  **DATA-SPECIFIC INFORMATION FOR: 6 files with names ERA5__**_100mm_TROTS4.txt**

1. Description: ERA5 climate data for all 483 sites. Used as input for all SEA analyses. File names indicate climate variable (pre=P, vpd=VPD, cwd=CWD).
2. Number of variables: 483
3. Number of rows: 81
4. Variable List:

* Each column contains the seasonal climate data for one chronology (site). Column name corresponds to the SiteID from metadata_FINAL.csv
* Each row corresponds to a year, with the row name representing the calendar year

  **DATA-SPECIFIC INFORMATION FOR: 6 files with names TerraClimate__**_100mm_TROTS4.txt**

1. Description: TerraClimate climate data for all 483 sites. Used as input for all SEA analyses. File names indicate climate variable (pre=P, vpd=VPD, cwd=CWD).
2. Number of variables: 483
3. Number of rows: 63
4. Variable List:

* Each column contains the seasonal climate data for one chronology (site). Column name corresponds to the SiteID from metadata_FINAL.csv
* Each row corresponds to a year, with the row name representing the calendar year

  DATA-SPECIFIC INFORMATION FOR: 4 files with names Tropics_detrended_sitechronologies_biweightmean_trots2_**spline.txt

1. Description: Mean tree-ring width chronologies for all sites, from 1930, for 4 detrending options. File names contain the number of years set for the 50% cutoff of the spline used for detrending (20, 30, 40, 50).
2. Number of variables: 484
3. Number of rows: 91
4. Variable List:

* Each column contains the RWI values for one chronology (site). The column name corresponds to the SiteID from metadata_FINAL.csv
* Each row corresponds to a calendar year, with the row name representing the calendar year
  DATA-SPECIFIC INFORMATION FOR: 12 files with names
  Interp_***lag0_Wet seasonseason_climspaceannual_Gymnos_final.tif
  and 12 with names
  Interp***_lag0_Wet seasonseason_climspaceannual_Gymnos_SD_final.tif

1. Description: First set of 12 files: SpatRasters with mean climatic interpolation models for all combinations of drought type (P, VPD, CWD), season (dry, wet), and taxonomic group (Angiosperms, Gymnosperms). Second set of 12 files containing SD in the file names: SpatRasters of SD from interpolation analyses, for the same combinations. Used to produce Figure 4 and SI Figures S14, S15, and S17.
2. Number of variables: 1
3. Number of rows, columns: 400, 310
4. Variable List:

* Either mean or SD of interpolated RWI anomaly
  DATA-SPECIFIC INFORMATION FOR: 2 files with names
  Map_P_lag0_Dryseasonseason_geospace_***elev10perc_juniperrangeAfrica_final.tif
  and 2 files with names
  Map_P_lag0_Dryseasonseason_geospace***_elev10perc_juniperrangeAfrica_SD_final.tif

1. Description: First set of 2 files: SpatRasters with maps of mean growth anomalies from interpolation models for dry-season P droughts for Gymnosperms and wet-season P droughts for Angiosperms. Second set of 2 files: same, but now containing SD values. Used to produce Figure 4 and SI Figure S17.
2. Number of variables: 1
3. Number of rows, columns: 720, 4320
4. Variable List:

* Either mean or SD of interpolated RWI anomaly
  DATA-SPECIFIC INFORMATION FOR: 483 netcdf files with names consisting of SiteID followed by \_detrended.nc

1. Description: Each file contains all detrended ring-width series belonging to one chronology, using 4 detrending methods (20, 30, 40, and 50-year splines). Used to produce series-level SEA analyses, shown in SI Figure S8. Dimensions in netcdf represent: years, series, splines
2. Number of variables: 1
3. Number of dimensions: first dimension varies across chronologies depending on the duration of the chronology, second dimension varies across chronologies depending on the number of series included, third dimension equals 4, representing splines of 20, 30, 40, and 50 years.
4. Variable List:

* RWI of individual series
  DATA-SPECIFIC INFORMATION FOR: 483 netcdf files with names consisting of SiteID followed by \_detrended.nc

1. Description: Each file contains all detrended ring-width series belonging to one chronology, using 4 detrending methods (20, 30, 40, and 50-year splines). Used to produce series-level SEA analyses, shown in SI Figure S8. Dimensions in netcdf represent: years, series, splines
2. Number of variables: 1
3. Number of dimensions: first dimension varies across chronologies depending on the duration of the chronology, second dimension varies across chronologies depending on the number of series included, third dimension equals 4, representing splines of 20, 30, 40, and 50 years.
4. Variable List:

* RWI of individual series
  CODE OVERVIEW

1. File List:
   A) 1_Chronologies_TropicalTreeRingNetwork.R
   B) 2_Climate_TropicalTreeRingNetwork.R
   C) 3_SEA_TropicalTreeRingNetwork.R
   D) 4_Interpolation_TropicalTreeRingNetwork.R
   E) 5_Figures_TropicalTreeRingNetwork.R
   F) 6_SIFigures_TropicalTreeRingNetwork.R
2. Software to run code: R
   CODE-SPECIFIC INFORMATION FOR: 1_Chronologies_TropicalTreeRingNetwork.R
3. Description: Builds chronologies from raw ring-width data, applying different detrending methods.
4. Uses: Raw ring-width data available at the ITRDB through links in metadata_FINAL.csv
5. Produces: 5 txt files with chronologies, differing in detrending method
   CODE-SPECIFIC INFORMATION FOR: 2_Climate_TropicalTreeRingNetwork.R
6. Description: Collects climate data, calculates derived climate variables (CWD, VPD), and calculates seasonal values (dry and wet season).
7. Uses: Raw climate data from three climate products: CRU, TerraClimate, and ERA5.
8. Produces: 6 txt files per climate product, for 3 climate variables (P, VPD, CWD) and two seasons (dry, wet)
   CODE-SPECIFIC INFORMATION FOR: 3_SEA_TropicalTreeRingNetwork.R
9. Description: Conducts superposed epoch analysis for all chronologies, using various drought definitions (10%, 5%, 2-year droughts), at chronology and series level, and for three climate products (CRU, ERA5, TerraClimate).
10. Uses: Seasonal climate data, chronologies.
11. Produces: sea_1_standard.Rdata, sea_2_wet.csv, sea_3_series.csv, sea_4_splines.csv, sea_5_CRU_ERA.Rdata, sea_6_CRU_TC.Rdata, sea_7_fivepercent.csv, sea_8_repeated.csv
    CODE-SPECIFIC INFORMATION FOR: 4_Interpolation_TropicalTreeRingNetwork.R
12. Description: Interpolates growth anomalies in climate space and maps these interpolation results in geographic space.
13. Uses: RWI anomalies from sea_1_standard.Rdata
14. Produces: TIFs of mean & SD growth anomalies in climate & geographic space
    CODE-SPECIFIC INFORMATION FOR: 5_Figures_TropicalTreeRingNetwork.R
15. Description: Code to conduct statistical analyses and produce main Figures.
16. Uses: metadata_FINAL.csv, sea_1_standard.Rdata, sea_2_wet.Rdata, MODIS_tree_cover.tif
17. Produces: 4 main text Figures
    CODE-SPECIFIC INFORMATION FOR: 6_SIFigures_TropicalTreeRingNetwork.R
18. Description: Code to conduct statistical analyses and produce SI Figures, Tables, and Dataset S1.
19. Uses: metadata_FINAL.csv,sea_1_standard.Rdata, sea_2_wet.csv, sea_3_series.csv, sea_4_splines.csv, sea_5_CRU_ERA.Rdata, sea_6_CRU_TC.Rdata, sea_7_fivepercent.csv, sea_8_repeated.csv, MODIS_tree_cover.tif, S_mean_raster.tifmetadata_FINAL.csv, sea_1_standard.Rdata, sea_2_wet.Rdata, MODIS_tree_cover.tif, all interpolation GeoTiffs.
20. Produces: 18 SI Figures, 1 Table, 1 Dataset

