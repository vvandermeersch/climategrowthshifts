Title        : Productivity of riparian Populus forests: satellite assessment
               along a prairie river with an environmental flow regime

Author       : Oscar R. Zimmerman*, Stewart B. Rood*, Lawrence B. Flanagan*^
information    *Department of Biological Sciences, University of Lethbridge, 
                4401 University Drive W., Lethbridge, Alberta T1K 3M4, Canada 
               ^Corresponding author (E-mail: larry dot flanagan at uleth dot ca)

Abstract     : In semi-arid regions, the growth and survival of cottonwoods 
               (riparian Populus species) depend on river water supplementing the
               limited precipitation. Indicators of growth and productivity are 
               needed to assess how altered streamflow regimes on regulated rivers 
               impact cottonwood trees and the riparian forest ecosystems they 
               support. Satellite imagery from the Landsat program was used to
               make historical assessments of ecosystem productivity in a riparian 
               cottonwood forest along a regulated prairie river in southern 
               Alberta, Canada from 1984 to 2020, with an environmental flow 
               regime that increased the minimum flows implemented in 1993. A 
               version of the near-infrared reflectance of vegetation scaled with
               incoming sunlight (NIRvP) was calculated from Landsat images to 
               provide a proxy for primary production. NIRvP was validated 
               against gross primary production measurements from eddy covariance
               and cottonwood basal area increment measurements from tree ring 
               analyses. Streamflow and weather data were used to assess what 
               environmental conditions drive year-to-year variations in NIRvP.  

Associated   : This dataset compliments the following manuscript: 
publiations     
and datasets   Zimmerman, O. R., Rood, S. B., & Flanagan, L. B. (2022). 
               Productivity of riparian Populus forests: satellite assessment
               along a prairie river with an environmental flow regime. Ecosphere. 
	       
               hereby referred to as 'Zimmerman et al. (2022)'. 
               NOTE: Currently under review. This will be updated with a proper 
               citation upon acceptance/publication of the manuscript. 

               Some of the data in this dataset and the associated manuscript have
               been previously published. Details are provided below under
               'Methods'. 

Files        : zimmerman_2022_riparian_productivity_figure3-4.csv
               zimmerman_2022_riparian_productivity_figure5.csv
               zimmerman_2022_riparian_productivity_figure6.csv
               zimmerman_2022_riparian_productivity_figure7-8.csv  
               zimmerman_2022_riparian_productivity_figure9.csv
               zimmerman_2022_riparian_productivity_table1-2.csv
	       NOTE: Details provided in separate sections below. 

Location of  : 'Helen Schuler Nature Reserve' field site (49.702 degrees N, 
collected       112.863 degrees W), Lethbridge, Alberta, Canada   
data

Methods      : Images from the Landsat 5 Thematic Mapper (TM), Landsat 7 Enhanced
               Thematic Mapper (ETM+), and Landsat 8 Operational Land Imager (OLI)
               from 1984 to 2020 were combined. Colleciton 1 Level 1 data products
               were used, which are precision- and terrain-corrected to the highest
               standard and processed to surface reflectance. Images spanned two 
               World Reference System-2 (WRS-2) Paths/Rows (41/25 and 40/25). 
               The pixel quality assessment band was used to mask clouds and cloud 
               shadows, and snow pixels were flagged. Gains were applied that
               normalized TM and OLI reflectances to ETM+-like values. Images were 
               accessed from the Google Earth Engine 
               (https://earthengine.google.com/) and exported out of the Google 
               Cloud for further processing in MATLAB (Version 9.0, The Mathworks, 
               Natick, MA, USA; https://www.mathworks.com/products/matlab.html). 
               The near-infrared reflectance of vegetation was calculated from 
               spectral reflectances in the red and near-infrared bands. Temporal 
               changes in background reflectance (i.e. when leaf area index
               equalled zero) were accounted for in the NIRv calculation, using a
               minimum four-year lagged moving window. A simplified NIRvP was 
               calculated as the product of NIRv and daily potential incoming 
               sunlight (photosynthetically active radiation, PAR, in molar units).
               Spike-filtering and gap-filling were applied before fitting each
               pixel with a cubic spline (csaps function in MATLAB, setting the 
               smoothing parameter to 0.001) to retreive continuous daily NIRvP 
               values. To assess year-to-year patterns of NIRvP, two aggregated
               metrics were calculated. The growing season NIRvP was calculated as
               the daily sum from May to mid-October, with some exceptions to 
               these bound (e.g. years with snowfall during late spring or early 
               autumn). Late summer NIRvP was calculated as the average NIRvP from 
               mid-July to August. NIRvP values were spatially averaged over pixels
               contained in the studied cottonwood forest. This assessed area 
               contained a mixture of cottonwood canopy and understory herbaceous 
               vegetation. See Flanagan et al. (2017, 2019) and Rood et al. (2013) 
               for further details on the study site, and Zimmerman et al. (2022) 
               for further details on the NIRvP calculations. 

	       Gross primary production (GPP) was derived from measurements of net 
               ecosystem carbon dioxide exchange made using the eddy covariance 
               (EC) technique during the 2014, 2015, 2017, and 2018 growing 
               seasons (approximately May--September). The EC system consisted of
               a 3D sonic anemometer (SAT; CSAT3, Campbell Scientific) and a fast-
               response, enclosed, infrared gas analyzer (IRGA; LI7200, LI-COR) 
               mounted 22.5 m above the ground surface on a telescoping aluminium 
               tower. High frequency EC data were processed using EddyPro
               (Version 5.2.1, LI-COR;
               https://www.licor.com/env/support/EddyPro/software.html) and half-
               hourly data were further processed using custom MATLAB scripts. 
               See Flanagan et al. (2017) for further details on the 
               instrumentation, ancillary measurements, and data processing. 

               Basal area increments (BAI) were derived from tree rings collected
               from male and female narrowleaf cottonwoods (Populus angustifolia) 
               during April and May 2009 (n=31, combined). Trees were in healthy 
               condition and of relatively uniform size. BAI were calculated from 
               tree diameters at breast height and annual radial increments
               assuming a cylindrical trunk. Radial increments were measured from 
               three cores per tree, collected from N-, SE-, and SW- facing 
               sectors, and averaged for BAI calculations. BAI were averaged over
               all trees for a given year for the time period 1992--2008. See 
               Rood et al. (2013) for further details on the trees and
               measurements. 
               
	       Daily river dishcarges (streamflow) were accessed from the Water
               Survey of Canada website (https://wateroffice.ec.gc.ca/) for the
               Oldman River near Lethbridge (station 05AA033). Daily precipitaition
               and air temperature (minimums and maximums), aggregated to monthly 
               values, were obtained from the Daymet Version 4 model output, 
               accessed from the Google Earth Engine 
               (https://earthengine.google.com/). Air temperature and precipition 
               data were used to calculate a soil moisture index (SMI) following 
               Hogg et al. (2013). SMI values were scaled [0,1], which differed 
               slightly from its original usage by Hogg and others. 

               GPP data for 2014, 2015, and 2017 have been previously published in 
               Flanagan et al. (2017) and Yang et al. (2019). BAI data have been 
               previously published in Rood et al. (2013). 

	       REFERENCES: 

	       Flanagan, L. B., Orchard, T. E., Logie, G. S., Coburn, C. A., & 
                  Rood, S. B. (2017). Water use in a riparian cottonwood ecosystem:
                  Eddy covariance measurements and scaling along a river corridor. 
                  Agricultural and Forest Meteorology, 232, 332–348.
               Flanagan, L. B., Orchard, T. E., Tremel, T. N., & Rood, S. B.  
                  (2019). Using stable isotopes to quantify water sources for trees
                  and shrubs in a riparian cottonwood ecosystem in flood and 
                  drought years. Hydrological Processes, 33(24), 3070–3083.
               Hogg, E. H., Barr, A. G., & Black, T. A. (2013). A simple soil 
                  moisture index for representing multi-year drought impacts on 
                  aspen productivity in the western Canadian interior. Agricultural
                  and Forest Meteorology, 178, 173–182.      
	       Rood, S. B., Ball, D. J., Gill, K. M., Kaluthota, S., 
                  Letts, M. G., & Pearce, D. W. (2013). Hydrologic linkages 
                  between a climate oscillation, river flows, growth, and wood 
                  ∆13C of male and female cottonwood trees. Plant, Cell & 
                  Environment, 36(5), 984–993.
	       Yang, H., Rood, S. B., & Flanagan, L. B. (2019). Controls on 
                  ecosystem water-use and water-use efficiency: Insights from a 
                  comparison between grassland and riparian forest in the northern
                  Great Plains. Agricultural and Forest Meteorology, 271, 22–32.

Revision     : 2021-12-21 (yyyy-mm-dd): metadata (README.txt) and data (.csv's)
history        files created by Oscar Zimmerman.
               
-----------------------------------------------------------------------------------

File         : zimmerman_2022_riparian_productivity_figure3-4.csv 

Associated   : This data can be used to reproduce Figures 3 and 4 in 
figures,       Zimmerman et al. (2022) and the following associated statistical 
tables,        analyses: linear regression, correspondence analysis, analysis 
statistics,    of covariance. 
etc. 

Column       : Year                         Year
descriptions   Jday                         Julian day
               NIRvP_spline_mol_m-2_d-1     NIRvP daily cubic spline fit 
               NIRvP_spline_7d_mol_m-2_d-1  NIRvP 7-day average cubic spline fit 
               NIRv_spline                  NIRv daily cubic spline-fit 
               NIRv_spline_7d               NIRv 7-day average cubic spline fit 
               GPP_g_C_m-2_d-1              Gross primary production
               GPP_flag                     Flagged days (see notes)
               GPP_7d_g_C_m-2_d-1           Gross primary production 7-day average

Column       : Year                         Years
units          Jday                         Days
               NIRvP_spline_mol_m-2_d-1     Moles of PAR per meter squared per day
               NIRvP_spline_7d_mol_m-2_d-1  Moles of PAR per meter squared per day
               NIRv_spline                  Unitless
               NIRv_spline_7d               Unitless 
               GPP_g_C_m-2_d-1              Grams of carbon per meter squared per day
               GPP_flag                     Binary
               GPP_7d_g_C_m-2_d-1           Grams of carbon per meter squared per day

Notes        : Flagged days correspond to days with GPP data removed from 
               subsequent analyses (see Zimmerman et al. 2022). 1's corespond to 
               flags (0's are not flagged). All blank data cells set to 'NaN'.

-----------------------------------------------------------------------------------

File         : zimmerman_2022_riparian_productivity_figure5.csv 

Associated   : This data can be used to reproduce Figure 5 in Zimmerman et al.
figures,       (2022) and following associated statistical analysis: 
tables,        correspondence analysis. 
statistics,    
etc.  

Column       : Year                        Years
descriptions   BAI_cm2                     Basal area increment
               NIRvP_GS_mol_m-2            Growing season integrated NIRvP 
               NIRv_GS_d                   Growing season integrated NIRv 

Column       : Year                        Years
units          BAI_cm2                     Centimeters squared 
               NIRvP_GS_mol_m-2            Moles of PAR per meter squared
               NIRv_GS_d                   Days 

-----------------------------------------------------------------------------------

File         : zimmerman_2022_riparian_productivity_figure6.csv 

Associated   : This data can be used to reproduce Figure 6 in
figures,       Zimmerman et al. (2022).
tables,         
statistics,    
etc.      

Column       : Year                  Year
descriptions   Q_GS_m3_s-1           Streamflow GS average
               Q_GS3_m3_s-1          Streamflow GS 3-year running average 
               P_GS_mm               Precipitation GS 
               P_GS3_mm              Precipitation GS 3-year running average 
               T_GS_degC             Air temperature GS average
               T_GS3_degC            Air temperature GS 3-year running average
               SMI_GS                Soil moisture index GS average 
               SMI_GS3               Soil moisture index 3-year running average 

Column       : Year                  Years
units          Q_GS_m3_s-1           Meters cubed per second
               Q_GS_m3_s-1           Meters cubed per second
               P_GS_mm               Millimeters
               P_GS3_mm              Millimeters
               T_GS_degC             Degrees Celsius 
               T_GS3_degC            Degrees Celsius 
               SMI_GS                Unitless 
               SMI_GS3               Unitless 

Notes        : GS = growing season (May--September). Running averages are centered.

-----------------------------------------------------------------------------------

File         : zimmerman_2022_riparian_productivity_figure7-8.csv 

Associated   : This data can be used to reproduce Figures 7 and 8 in 
figures,       Zimmerman et al. (2022).
tables,         
statistics,    
etc. 

Column       : Year                        Year  
descriptions   Jday                        Julian day 
               NIRvP_mol_m-2_d-1           NIRvP spatial average over all pixels
               NIRvP_spline_mol_m-2_d-1    NIRvP daily cubic spline fit 
               Q_m3_s-1                    Streamflow daily average

Column       : Year                        Years
units          Jday                        Days
               NIRvP_mol_m-2_d-1           Moles of PAR per meter squared per day  
               NIRvP_spline_mol_m-2_d-1    Moles of PAR per meter squared per day 
               Q_m3_s-1                    Meters cubed per second

Notes        : Note that cubic spline fits of NIRvP were calculated as the average 
               of splines fit to each pixel, not as the spline fit to the average 
               of all pixels. So NIRvP (spatial average over all pixels, col3) is  
               only used for visualization purposes to give a sense of available
               imagery. All blank data cells set to 'NaN'.

-----------------------------------------------------------------------------------

File         : zimmerman_2022_riparian_productivity_figure9.csv 

Associated   : This data can be used to reproduce Figure 9 in 
figures,       Zimmerman et al. (2022).
tables,         
statistics,    
etc. 

Column       : Year                   Year
descriptions   NIRvP_GS_mol_m-2       Growing season integrated NIRvP 
               NIRvP_LS_mol_m-2_d-1   Late summer NIRvP 
               Qy2_GS_m3_s-1          Streamflow 2-year GS average
               Qy2_annual_m3_s-1      Streamflow 2-year annual average 
               SMI_GS                 Soil moisture index GS average
               SMI_annual             Soil moisture index annual average

Column       : Year                   Years
units          NIRvP_GS_mol_m-2       Moles of PAR per meter squared
               NIRvP_LS_mol_m-2_d-1   Moles of PAR per meter squared per day 
               Qy2_GS_m3_s-1          Meters cubed per second
               Qy2_annual_m3_s-1      Meters cubed per second
               SMI_GS                 Unitless
               SMI_annual             Unitless 
     
Notes        : GS = growing season (May--September). Annual = prior October to 
               current September. 2-year averages calculated from the current year
               and prior year weighted one-half. 

-----------------------------------------------------------------------------------

File detail  : zimmerman_2022_riparian_productivity_table1-2.csv 

Associated   : This data can be used to reproduce Tables 1 and 2 in 
figures,       Zimmerman et al. (2022) which display outputs of the following 
tables,        statistical analyses: correlation analysis, regression tree  
statistics,    analysis. 
etc. 

Column       : Year                  Year
descriptions : Qy_ES_m3_s-1          Streamflow ES average (current year)
               Qy_LS_m3_s-1          Streamflow LS average (current year)
               Qy_GS_m3_s-1          Streamflow GS average (current year)
               Qy2_GS_m3_s-1         Streamflow 2-year GS average (current year) 
               Qy_annual_m3_s-1      Streamflow annual average (current year)
	       Qy2_annual_m3_s-1     Streamflow 2-year annual average (current year)
               Qyp_ES_m3_s-1         Streamflow ES avearge (prior year)
               Qyp_LS_m3_s-1         Streamflow LS average (prior year)
               Qyp_GS_m3_s-1         Streamflow GS average (prior year)
               Qyp_annual_m3_s-1     Streamflow annual average (prior year)
               Py_ES_mm              Precipitation ES (current year)
               Py_LS_mm              Precipitation LS (current year)
               Py_GS_mm              Precipitation GS (current year)
               Py_annual_mm          Precipitation annual (current year)
               Pyp_ES_mm             Precipitation ES (prior year)
               Pyp_LS_mm             Precipitation LS (prior year)
               Pyp_GS_mm             Precipitation GS (prior year)
               Pyp_annual_mm         Precipitation annual (prior year)
               Ty_ES_degC            Air temperature ES average (current year)
               Ty_LS_degC            Air temperature LS average (current year)
               Ty_GS_degC            Air temperature GS average (current year)
               Ty_annual_degC        Air temperature annual average (current year)
               Typ_ES_degC           Air temperature ES average (prior year)
               Typ_LS_degC           Air temperature LS average (prior year)
               Typ_GS_degC           Air temperature GS average (prior year)
               Typ_annual_degC       Air temperature annual average (prior year)
               SMI_ES                Soil moisture index ES average (current year)
               SMI_LS                Soil moisture index LS average (current year)
               SMI_GS                Soil moisture index GS average (current year)
               SMI_annual            Soil moisture index annual average (current year)
               SMIyp_ES              Soil moisture index ES average (prior year)
               SMIyp_LS              Soil moisture index LS average (prior year)
               SMIyp_GS              Soil moisture index GS average (prior year)
               SMIyp_annual          Soil moisture index annual average (prior year)   
               EnvFlows              Environmental flows 
               NIRvP_GS_mol_m-2      Growing season integrated NIRvP
               NIRvP_LS_mol_m-2_d-1  Late summer NIRvP 

Column       : Year                  Years
units          Columns 2--11         Meters cubed per second
               Columns 12--19        Millimeters 
               Columns 20--27        Degrees Celsius
               Columns 28--35        Unitless
               EnvFlows              Binary 
               NIRvP_GS_mol_m-2      Moles of PAR per meter squared
               NIRvP_LS_mol_m-2_d-1  Moles of PAR per meter squared per second

Notes        : ES = early season (May--June). LS = late season (July--August).
               GS = growing season (May--September). Annual = prior October to 
               current September. 2-year averages calculated from the current year
               and prior year weighted one-half. 