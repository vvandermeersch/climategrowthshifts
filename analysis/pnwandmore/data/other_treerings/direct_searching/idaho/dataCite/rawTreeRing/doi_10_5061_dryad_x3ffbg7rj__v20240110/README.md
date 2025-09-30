# Data from: A species' response to spatial climatic variation does not predict its response to climate change

<https://doi.org/10.5061/dryad.x3ffbg7rj>

This dataset contains tree ring timeseries and associated historical climate data used in the manuscript titled "A species' response to spatial climatic variation does not predict its response to climate change", by Daniel L. Perret, Margaret E. K. Evans, and Dov F. Sax. See Methods for description of data collection procedures. In brief, representative samples were collected from populations of ponderosa pine in the western United States that span the breadth of the climate space occupied by the species.

## Description of the data and file structure

The data are contained in a single .csv file. This file is a dataframe with 36,608 rows and 32 columns. Each row corresponds to an individual annual growth increment for some tree during some year in some location. The column labels and their meanings are thus:

*Series* - unique identifier for the tree core sample the data came from. Consists of two parts *site_tree* (e.g., "a3r1_10" indicates the data come from tree 10 sampled in site a3r1).
*Site, Tree* - information contained in *Series*
*Year* - the crossdated year assigned to the measurement; the year the growth ring was formed in
*RW* - the width (in mm) of the growth ring
*dch_mm* - the diameter at coring height of the tree (in mm) as measured in 2018 at sample collection
*Tree.Size* - the back-calculated diameter of the tree (in mm) in year *Year-1*
*Basal.Area* - *Tree.Size* converted to area (mm^2)
*next.Basal.Area* - *Basal.Area* in the subsequent *Year*
*BAI* - basal area increment, the difference between next.Basal.Area and Basal.Area
*x, y* - WGS84 latitude and longitude coordinates for the plot center associated with the Tree
*bio01* - mean annual temperature associated with the year and location indicated (PRISM)
*bio12* - cumulative annual precipitation associated with the year and location indicated (PRISM)
*bio01.norm* - long-term mean annual temperature (over timeseries length) associated with location indicated (PRISM)
*bio12.norm* - long-term mean annual precipitation (over timeseries length) associated with location indicated (PRISM)
*p.jj.tmax* - *so.tmax* - monthly maximum temperature summarized in two to four-month seasonal periods, indicated by the letters preceding "."; for example, "so_tmax" indicates the maximum temperature in september and october during the indicated Year in the indicated location. Variables that begin with "p." correspond with the maximum temperature in the indicated seasonal period in the previous year. 'NA' values indicate that this variable could not be calculated for the given year because it required data from before the timeseries began.
*p.jj.ppt* - *so.ppt* - monthly maximum temperature summarized in two to four-month seasonal periods, indicated by the letters preceding "."; for example, "so_tmax" indicates the maximum temperature in september and october during the indicated Year in the indicated location. Variables that begin with "p." correspond with the maximum temperature in the indicated seasonal period in the previous year. 'NA' values indicate that this variable could not be calculated for the given year because it required data from before the timeseries began.

## Sharing/Access information

These data and associated summaries are also available at the corresponding author's GitHub repository: <https://github.com/daniel-perret/PIPO_treerings>

Climate data was derived from the following source:

PRISM ( prism.oregonstate.edu )

## Code/Software

All code necessary to manipulate/generate data and perform analyses described in the associated manuscript are available at the corresponding author's GitHub respository: <https://github.com/daniel-perret/PIPO_treerings>
