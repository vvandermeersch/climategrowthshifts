############
Description of climate files in treerings_ailene

############
See also 2011 Ettinger, Ford and HRL (Ecology) paper (and supplement)


####Climate data
Nine climate variables were used to explore climate - tree growth relationships. These variables were calculated for each of the nine stands from the Longmire climate record (daily values of temperature and precipitation) using gridded PRISM climate data to estimate offsets (differences in temperature and precipitation between Longmire and stand locations). These estimated stand specific daily values were then used to calculate yearly average temperature and cumulative precipitation (MAT, PTT) as well as growing season and dormant season temperature and precipitation (GST, DST, GPT, DPT) for each water year (October 1 - September 30). Daily temperature values were also used to calculate growing degree days (GDD), and daily temperature and precipitation values were run through the Anderson snow model to estimate the maximum size of the snowpack (SWE) and snow duration (number of days with snow covered ground). 

1. tree_plot_climate_temp.csv - contains monthly (average) temperature values (in Celsius) from 1909 to 2009 (rows), for all 15 PSP stands (rows). These data were used to calculate annual and seasonal values.

2. tree_plot_climate_precip.csv - contains monthly precipitation (in cm) from 1909 to 2009 in columns, for all PSP stands (rows). These data were used to calculate annual and seasonal values.

3. tree_plot_climate_grdd5.csv - contains growing degree day data (the annual sum of daily mean temperatures for days with mean temperatures above 5C) for each stand (column) and each year (row) from 1914-2009.

4. tree_plot_climate_snowdur.csv - contains the duration of snowcover (in days) for each stand (column) and each year (row) from 1914-2009. 

5. tree_plot_climate_swe.csv - contains the maximum size of the snowpack (in mm snow water equivalent) for each stand (column) and each year (row) from 1914-2009. 






