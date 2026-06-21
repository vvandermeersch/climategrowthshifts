## Description of climate files in treerings_ailene
Started at BIRS in March 2025 by Janneke HRL with updates by Lizzie

See also 2011 Ettinger et al. 2011 (https://doi.org/10.1890/10-1639.1), Ford and HRL 2016 (https://cdnsciencepub.com/doi/10.1139/cjfr-2016-0188) papers


## Climate data
Nine climate variables were used to explore climate - tree growth relationships. These variables were calculated for each of the nine stands from the Longmire climate record (daily values of temperature and precipitation) using gridded PRISM climate data to estimate offsets (differences in temperature and precipitation between Longmire and stand locations). These estimated stand specific daily values were then used to calculate yearly average temperature and cumulative precipitation (MAT, PTT) as well as growing season and dormant season temperature and precipitation (GST, DST, GPT, DPT) for each water year (October 1 - September 30). Daily temperature values were also used to calculate growing degree days (GDD), and daily temperature and precipitation values were run through the Anderson snow model to estimate the maximum size of the snowpack (SWE) and snow duration (number of days with snow covered ground). 

1. tree_plot_climate_temp.csv - contains monthly (average) temperature values (in Celsius) from 1909 to 2009 (rows), for all 15 PSP stands (rows). These data were used to calculate annual and seasonal values.

2. tree_plot_climate_precip.csv - contains monthly precipitation (in cm) from 1909 to 2009 in columns, for all PSP stands (rows). These data were used to calculate annual and seasonal values.

3. tree_plot_climate_grdd5.csv - contains growing degree day data (the annual sum of daily mean temperatures for days with mean temperatures above 5C) for each stand (column) and each year (row) from 1914-2009.

4. tree_plot_climate_snowdur.csv - contains the duration of snowcover (in days) for each stand (column) and each year (row) from 1914-2009. 

5. tree_plot_climate_swe.csv - contains the maximum size of the snowpack (in mm snow water equivalent) for each stand (column) and each year (row) from 1914-2009. 


## Follow up with Ailene about neighborhood data June 2026

I originally wrote: 

I now write another project. With Victor (cc-ed) and Devina (an undergraduate in the lab, cc-ed) we have been using your attached data and PhD tree cores to try a new modeling approach to competition. We are currently deciding if we need to add year to the model. For that, we need to know a little more of how the data were collected and cannot find it in the papers. My understanding is that you cored outside the PSP plots and collected your own neighborhood data (attached), correct? If so, why are there multiple years of neighboorhood data? See:

https://github.com/vvandermeersch/climategrowthshifts/issues/123#issuecomment-4712348162

However, there is only ever ONE year for any ONE focal tree. 

And Ailene replied:

 Hello Lizzie, Victor and Devina-
So great that you're using these data! We did not collect all the neighbor data in a single year- it took us 3-4 field seasons to collect all the neighbor data for the cored trees as it is very time consuming. I believe we started by coring 10 trees per elevation on the south side in 2008, then expanded to the east side and to 20 trees per elevation (so coring 10 more at sites we had only cored 10). We decided to collect neighbor data starting in the second year of data collection, I believe, so went back to all previously cored trees that we had NOT collected neighbor data at. 
Hope that makes sense- let me know if it does not or if you have any follow up questions.
Ailene



