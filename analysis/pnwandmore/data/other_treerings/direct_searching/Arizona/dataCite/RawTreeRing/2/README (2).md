# Tree-ring width, and tree-ring stable carbon isotope data from SW Alaska

Dataset DOI: [10.5061/dryad.9w0vt4bsm](10.5061/dryad.9w0vt4bsm)

## Description of the data and file structure

The study sites comprise a subset of a larger network of monitoring plots on the northern Alaska Peninsula, on and adjacent to National Park Service lands (Csank et al. 2016; Sherriff et al. 2017). At each site, we preferentially cored a minimum of 20 mature, large diameter canopy trees at a distance of at least 10 m from one another. At forest sites we sampled trees >20 cm DBH and >10 m in height. At woodland sites, where trees were smaller, we sampled trees >8 cm DBH. Additional details regarding field methods can be found in previously published studies from these sites (Csank et al. 2016; and Sherriff et al. 2017). 

Cores were collected at 30-50 cm above the root crown and crossdated using standard dendrochronological techniques (Stokes and Smiley 1996). Tree-ring widths for each core were counted using a binocular scope, measured to 0.001 mm using WinDendro (Regent Instruments, Inc.) at a scan resolution of 600-1200 dpi, and were visually and statistically crossdated using the COFECHA program (Holmes et al. 1983). We generally found very consistent inter-series growth patterns within sites, with inter-series correlations averaging 0.47 across sites. 

Prior to detrending long-term growth trends and averaging the individual indices, we grouped the cores from each site according to their recent growth patterns to develop separate chronologies for the following: (1) trees that exhibited at least a 2-fold increase in growth during the most recent 10 years, relative to growth in the prior decade; and (2) trees that showed reduced or stationary growth over the same period (n = 6-10 cores per group). Two sites affected by beetle kill included separate chronologies from live and dead trees, both of which were identified as having reduced or stationary growth. A third beetle kill site included cores from dead trees only. The groupings allowed us to differentiate between trees with recent growth increases well above long-term growth trends (i.e., >2-fold increase in the most recent decade) from other trees at the same site that showed declining or stationary growth. It is possible that some trees identified as having stationary growth could have been grouped with trees showing increasing growth if we had used a lower threshold (e.g. 1.5x), but we believe that the 2x threshold effectively captures the variation in growth patterns that we observed at our sites. Five sites had chronologies that showed only increasing growth, two sites had only stationary growth, and four sites included separate chronologies that showed both increasing and stationary growth. 

For more details please see Miller et al. (2025) Ecology

### Files and variables

#### File: Miller\_et\_al\_2025\_Ecology\_AK\_Spruce\_Data.xlsx

**Description:** Excel worksheet contains 4 sheets. The first is metadata including site locations, species, coordinates and time period covered. The second worksheet is RWI for each site, the third is the raw carbon isotope data and the fourth is the corrected carbon isotope discrimination data. Data were corrected using the Point Barrow d13CO2 and PCO2 data. Missing data is indicated as 'NaN'

##### Variables

* RWI (dimentionless)
* d13C (raw) (per mil)
* D13C (discrimination) (per mil)

## Code/software

The R statistical program was used with the packages dplR. 

## Access information

Other publicly accessible locations of the data:

* These data will also be uploaded to the ITRDB

